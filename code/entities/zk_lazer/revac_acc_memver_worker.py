#!/usr/bin/env python3
"""Persistent binary local bridge for AccGenpp MemVer Labrador proofs.

The C protocol sends a short dynamic request keyed by a 32-byte MemVer public
parameter identifier.  On a cache miss the worker asks for the binary public
parameter blob once, then retains only that *public* data in an LRU cache.
Private UAV witnesses are decoded for one prove request and discarded.

This is intentionally a same-device Unix-domain-socket service.  It is not a
remote proving service and must never be used to send a UAV witness to TA/GCS.
"""

from __future__ import annotations

import argparse
import contextlib
import gc
import hashlib
import io
import os
import select
import signal
import socket
import struct
import sys
import time
import zlib
from collections import OrderedDict
from typing import Any

from revac_acc_memver_import_labrador import (
    LiveLabrador38Params,
    MemverStaticTemplate,
    build_labrador,
    labrador38_proof_is_compact,
    labrador38_proof_params_bytes,
    labrador38_proof_profile_id,
    memver_labrador_profile_id,
    verify_labrador_live_statement,
)
from revac_profile_generated import (
    D,
    GADGET_BASE,
    GADGET_DIM,
    GADGET_SIGNED,
    HASH_BYTES,
    HANDLE_BITS,
    N,
    LABRADOR_WORKER_PROTOCOL_VERSION,
    Q,
    SEED_BYTES,
    TRUNC_K_L,
)


MAGIC = b"RVMVWB" + f"{LABRADOR_WORKER_PROTOCOL_VERSION:02d}".encode("ascii")
MAX_FRAME = 64 * 1024 * 1024
REFRESH_BITS_MAX = 192
REFRESH_BYTES_MAX = REFRESH_BITS_MAX // 8
PARAM_CACHE_LIMIT = 4
DOMAIN_SEPARATOR_DS = 5


class CompactRefreshMatrixBank:
    """Contiguous int32 storage for A_0,...,A_l.

    Keeping 6.3 million coefficients as Python ``int`` objects was the main
    reason the persistent worker grew past 2.7 GiB.  This wrapper owns either
    the public-blob bytes or a native CFFI allocation and materializes one
    polynomial only while the fixed statement is compiled.
    """

    def __init__(self, coeffs: Any, owner: Any, matrices: int):
        self.coeffs = coeffs
        self.owner = owner
        self.matrices = matrices

    def poly(self, matrix: int, block: int, row: int, col: int) -> list[int]:
        if not (0 <= matrix < self.matrices and
                0 <= block < GADGET_DIM and 0 <= row < D and 0 <= col < D):
            raise IndexError("refresh matrix polynomial index out of range")
        start = ((((matrix * GADGET_DIM + block) * D + row) * D + col) * N)
        return list(self.coeffs[start:start + N])


class CompactPolyTensor:
    """Lazy tensor view over contiguous int32 polynomial coefficients."""

    def __init__(self, coeffs: Any, owner: Any, shape: tuple[int, ...],
                 poly_offset: int = 0):
        if not shape or any(dimension <= 0 for dimension in shape):
            raise ValueError("compact polynomial tensor has an invalid shape")
        self.coeffs = coeffs
        self.owner = owner
        self.shape = shape
        self.poly_offset = poly_offset

    def __len__(self) -> int:
        return self.shape[0]

    def __getitem__(self, index: int) -> Any:
        if index < 0:
            index += self.shape[0]
        if not 0 <= index < self.shape[0]:
            raise IndexError(index)
        stride = 1
        for dimension in self.shape[1:]:
            stride *= dimension
        offset = self.poly_offset + index * stride
        if len(self.shape) == 1:
            start = offset * N
            return list(self.coeffs[start:start + N])
        return CompactPolyTensor(
            self.coeffs, self.owner, self.shape[1:], offset)


class CompactI64Polys:
    """Zero-copy view of C-generated exact-lift polynomials."""

    def __init__(self, coeffs: Any, owner: bytes, count: int):
        self.coeffs = coeffs
        self.owner = owner
        self.count = count

    def __len__(self) -> int:
        return self.count

    def __getitem__(self, index: int) -> list[int]:
        if not 0 <= index < self.count:
            raise IndexError(index)
        start = index * N
        return list(self.coeffs[start:start + N])


def _uniform_poly(seed: bytes, row: int, col: int) -> list[int]:
    """Match poly_q_uniform on the little-endian x86 deployment target."""
    absorbed = seed + struct.pack("<IQQ", DOMAIN_SEPARATOR_DS, row, col)
    # 2048 bytes is well above the expected ~826 bytes for 256 accepted
    # 19-bit coefficients at q=405769.  Grow deterministically if needed.
    out_len = 2048
    while True:
        stream = hashlib.shake_128(absorbed).digest(out_len)
        coeffs: list[int] = []
        off = 0
        while off + 5 <= len(stream) and len(coeffs) < N:
            word = int.from_bytes(stream[off:off + 5], "little")
            lo = word & ((1 << 19) - 1)
            hi = (word >> 20) & ((1 << 19) - 1)
            if lo < Q:
                coeffs.append(lo)
            if len(coeffs) < N and hi < Q:
                coeffs.append(hi)
            off += 5
        if len(coeffs) == N:
            return coeffs
        out_len *= 2


def _expand_refresh_matrix_bank(seed: bytes, refresh_bits: int) -> CompactRefreshMatrixBank:
    """Expand A_0,...,A_l from the canonical 32-byte public seed."""
    if os.environ.get("REVAC_REFRESH_SEED_EXPAND", "native") != "python":
        from lazer import ffi, lib

        matrix_count = refresh_bits + 1
        poly_count = matrix_count * GADGET_DIM * D * D
        coeffs = ffi.new("int32_t[]", poly_count * N)
        seed_buf = ffi.new("uint8_t[32]", seed)
        error = lib.labrador38_py_seeded_uniform_matrix_bank(
            coeffs, seed_buf, matrix_count, GADGET_DIM, D, D, N,
            Q, 19, DOMAIN_SEPARATOR_DS)
        if error != 0:
            raise ValueError(f"native refresh matrix expansion failed: {error}")
        return CompactRefreshMatrixBank(coeffs, coeffs, matrix_count)

    # Portable reference path retained for cross-checks and platforms where
    # the native bridge primitive is unavailable.
    matrices: list[Any] = []
    for matrix_index in range(refresh_bits + 1):
        for block in range(GADGET_DIM):
            block_seed = hashlib.shake_256(
                seed + struct.pack("<QQ", matrix_index, block)
            ).digest(SEED_BYTES)
            matrices.append([
                [_uniform_poly(block_seed, row, col) for col in range(D)]
                for row in range(D)
            ])
    from lazer import ffi
    flat = [coefficient for matrix in matrices for row in matrix
            for poly in row for coefficient in poly]
    coeffs = ffi.new("int32_t[]", flat)
    return CompactRefreshMatrixBank(coeffs, coeffs, refresh_bits + 1)


def recv_exact(conn: socket.socket, nbytes: int) -> bytes:
    out = bytearray()
    while len(out) < nbytes:
        part = conn.recv(nbytes - len(out))
        if not part:
            raise ConnectionError("peer closed the local Labrador bridge")
        out.extend(part)
    return bytes(out)


def recv_frame(conn: socket.socket) -> bytes:
    length = struct.unpack("<I", recv_exact(conn, 4))[0]
    if length == 0 or length > MAX_FRAME:
        raise ValueError(f"invalid worker frame length: {length}")
    return recv_exact(conn, length)


def send_frame(conn: socket.socket, payload: bytes) -> None:
    if not payload or len(payload) > MAX_FRAME:
        raise ValueError("invalid worker response length")
    conn.sendall(struct.pack("<I", len(payload)) + payload)


def read_i32_polys(raw: memoryview, off: int, count: int) -> tuple[list[list[int]], int]:
    nbytes = count * N * 4
    if off + nbytes > len(raw):
        raise ValueError("truncated i32 polynomial data")
    values = struct.unpack_from(f"<{count * N}i", raw, off)
    polys = [list(values[i * N:(i + 1) * N]) for i in range(count)]
    return polys, off + nbytes


def read_i64_polys(raw: memoryview, off: int, count: int) -> tuple[CompactI64Polys, int]:
    nbytes = count * N * 8
    if off + nbytes > len(raw):
        raise ValueError("truncated i64 polynomial data")
    from lazer import ffi
    owner = bytes(raw[off:off + nbytes])
    coeffs = ffi.from_buffer("int64_t[]", owner)
    return CompactI64Polys(coeffs, owner, count), off + nbytes


def decode_public(raw: bytes) -> dict[str, Any]:
    data = memoryview(raw)
    if len(data) < len(MAGIC) + SEED_BYTES or bytes(data[:len(MAGIC)]) != MAGIC:
        raise ValueError("bad MemVer public-parameter binary magic")
    off = len(MAGIC)
    seed = bytes(data[off:off + SEED_BYTES])
    off += SEED_BYTES
    if off + 4 + SEED_BYTES + 1 > len(data):
        raise ValueError("truncated randomized-refresh public header")
    refresh_bits = struct.unpack_from("<I", data, off)[0]
    off += 4
    refresh_seed = bytes(data[off:off + SEED_BYTES])
    off += SEED_BYTES
    refresh_blob_mode = int(data[off])
    off += 1
    if refresh_bits not in (0, 128, 192):
        raise ValueError("unsupported randomized-refresh bit length")
    if refresh_blob_mode not in (0, 1):
        raise ValueError("unsupported randomized-refresh public blob mode")
    from lazer import ffi
    base_matrices = 1 + TRUNC_K_L + HANDLE_BITS * GADGET_DIM + 2
    base_poly_count = base_matrices * D * D
    base_nbytes = base_poly_count * N * 4
    if off + base_nbytes > len(data):
        raise ValueError("truncated MemVer base public matrix bank")
    base_owner = bytes(data[off:off + base_nbytes])
    base_coeffs = ffi.from_buffer("int32_t[]", base_owner)
    off += base_nbytes
    poly_cursor = 0

    def tensor(shape: tuple[int, ...]) -> CompactPolyTensor:
        nonlocal poly_cursor
        result = CompactPolyTensor(base_coeffs, base_owner, shape, poly_cursor)
        count = 1
        for dimension in shape:
            count *= dimension
        poly_cursor += count
        return result

    atrap = tensor((D, D))
    trunc = tensor((TRUNC_K_L, D, D))
    ba = tensor((HANDLE_BITS, GADGET_DIM, D, D))
    link_a = tensor((D, D))
    link_d = tensor((D, D))
    if poly_cursor != base_poly_count:
        raise ValueError("bad MemVer base public matrix shape")
    if refresh_bits == 0:
        if refresh_blob_mode != 1 or off != len(data):
            raise ValueError("disabled refresh must use an empty seeded bank")
        refresh_A = None
    elif refresh_blob_mode == 0:
        refresh_nbytes = ((refresh_bits + 1) * GADGET_DIM * D * D * N * 4)
        if off + refresh_nbytes != len(data):
            raise ValueError("bad explicit refresh matrix-bank length")
        owner = bytes(data[off:off + refresh_nbytes])
        coeffs = ffi.from_buffer("int32_t[]", owner)
        refresh_A = CompactRefreshMatrixBank(
            coeffs, owner, refresh_bits + 1)
        off += refresh_nbytes
    else:
        refresh_A = _expand_refresh_matrix_bank(refresh_seed, refresh_bits)
    if off != len(data):
        raise ValueError("trailing bytes in MemVer public-parameter blob")
    return {
        "params": {
            "n": N, "d": D, "q": Q, "handle_bits": HANDLE_BITS,
            "gadget_base": GADGET_BASE, "gadget_dim": GADGET_DIM,
            "gadget_signed": int(GADGET_SIGNED), "trunc_k_l": TRUNC_K_L,
            "deletion_count": 0,
            "refresh_bits": refresh_bits,
        },
        "a_seed": seed,
        "Atrap": atrap,
        "trunc_B": trunc,
        "Ba": ba,
        "link_A": link_a,
        "link_D": link_d,
        "refresh_seed": refresh_seed,
        "refresh_A": refresh_A,
    }


def decode_dynamic(raw: bytes, public: dict[str, Any], mode: bytes) -> tuple[dict[str, Any], bytes | None]:
    data = memoryview(raw)
    if len(data) < 4:
        raise ValueError("truncated MemVer dynamic payload")
    deletion_count = struct.unpack_from("<I", data, 0)[0]
    off = 4
    acc, off = read_i32_polys(data, off, D)
    link_commitment, off = read_i32_polys(data, off, D)
    out = dict(public)
    out["params"] = dict(public["params"], deletion_count=deletion_count)
    out["acc"] = acc
    out["link_commitment"] = link_commitment
    if mode == b"V":
        if off + 4 > len(data):
            raise ValueError("truncated serialized proof length")
        proof_len = struct.unpack_from("<I", data, off)[0]
        off += 4
        if proof_len == 0 or off + proof_len != len(data):
            raise ValueError("invalid serialized proof length")
        return out, bytes(data[off:off + proof_len])
    if mode != b"P" or off + 4 > len(data):
        raise ValueError("invalid MemVer prove payload")
    handle = struct.unpack_from("<I", data, off)[0]
    off += 4
    link_r0, off = read_i32_polys(data, off, D)
    link_r1, off = read_i32_polys(data, off, D)
    left, off = read_i32_polys(data, off, D)
    v12, off = read_i32_polys(data, off, D)
    v2_flat, off = read_i32_polys(data, off, TRUNC_K_L * D)
    tag, off = read_i32_polys(data, off, 1)
    right_flat, off = read_i32_polys(data, off, HANDLE_BITS * GADGET_DIM * D)
    if off + 4 + REFRESH_BYTES_MAX > len(data):
        raise ValueError("truncated randomized-refresh witness header")
    refresh_bits = struct.unpack_from("<I", data, off)[0]
    off += 4
    refresh_r_raw = bytes(data[off:off + REFRESH_BYTES_MAX])
    off += REFRESH_BYTES_MAX
    if refresh_bits != int(public["params"]["refresh_bits"]):
        raise ValueError("randomized-refresh bit length mismatch")
    refresh_z_flat, off = read_i32_polys(data, off, GADGET_DIM * D)
    lift, off = read_i64_polys(data, off, D)
    link_lift, off = read_i64_polys(data, off, D)
    if off != len(data):
        raise ValueError("trailing bytes in MemVer prove payload")
    right = []
    pos = 0
    for _handle in range(HANDLE_BITS):
        gadgets = []
        for _gadget in range(GADGET_DIM):
            gadgets.append(right_flat[pos:pos + D])
            pos += D
        right.append(gadgets)
    out.update({
        "x": handle,
        "x_bits": [(handle >> bit) & 1 for bit in range(HANDLE_BITS)],
        "refresh_r_bits": [
            (refresh_r_raw[bit >> 3] >> (bit & 7)) & 1
            for bit in range(refresh_bits)
        ],
        "lift": lift,
        "link_lift": link_lift,
        "witness": {
            "left": left,
            "v12": v12,
            "v2": [v2_flat[i * D:(i + 1) * D] for i in range(TRUNC_K_L)],
            "tag": tag[0],
            "link_r": [link_r0, link_r1],
            "right": right,
            "refresh_z": [
                refresh_z_flat[j * D:(j + 1) * D]
                for j in range(GADGET_DIM)
            ],
        },
    })
    return out, None


def ensure_template(entry: dict[str, Any]) -> MemverStaticTemplate:
    template = entry.get("template")
    if template is None:
        template = MemverStaticTemplate(entry["public"])
        entry["template"] = template
    return template


def parameter_cache_key(data: dict[str, Any]) -> tuple[int, int, int, bool, str]:
    params = data["params"]
    return (int(params["deletion_count"]), int(params["gadget_base"]),
            int(params["gadget_dim"]), True, "grouped")


def close_entry(entry: dict[str, Any]) -> None:
    """Release every native/public object owned by one cached profile."""
    for params in entry.get("params", {}).values():
        params.close()
    entry.get("params", {}).clear()
    template = entry.get("template")
    if template is not None:
        # proof_statement owns its native statement and witness through its
        # destructor.  Removing the final template references lets CPython
        # release those objects before malloc_trim below.
        template.live_statements.clear()
        template.compiled_relation = None
        template.compiled_raw_relation = None
    entry["template"] = None
    entry["active_statement_profile"] = None
    entry.get("public", {}).clear()


def reset_global_comkey() -> None:
    """Release LaBRADOR's process-global commitment key at a safe point."""
    from lazer import lib
    lib.labrador38_py_comkey_reset()


def release_parameter_lease(entry: dict[str, Any]) -> None:
    """Optional low-memory policy for the large native Pack parameter tree."""
    if os.environ.get("REVAC_LABRADOR_PARAM_POLICY", "resident") != "ephemeral":
        return
    for params in entry.get("params", {}).values():
        params.close()
    entry.get("params", {}).clear()
    # Pack parameters reference the process-global commitment key.  In the
    # low-memory policy no live parameter object remains, so the key can also
    # be deterministically regenerated on the next request.
    from lazer import lib
    lib.labrador38_py_comkey_reset()


def compact_public_after_freeze(entry: dict[str, Any]) -> None:
    """Drop public Python matrix graphs once the native statement owns them."""
    template = entry.get("template")
    if template is None or not template.sealed:
        return
    public = entry["public"]
    for name in ("Atrap", "trunc_B", "Ba", "link_A", "link_D"):
        public.pop(name, None)


def trim_worker_heap() -> None:
    """Return compiler/proof temporaries to the OS at a policy safe point."""
    started = time.perf_counter()
    gc.collect()
    from lazer import lib
    lib.labrador38_py_malloc_trim()
    if os.environ.get("REVAC_LABRADOR_PROFILE", "0") != "0":
        print(
            "[Labrador worker 计时] 堆回收："
            f"{time.perf_counter() - started:.6f} 秒",
            file=sys.stderr,
            flush=True,
        )


def report_worker_memory(label: str) -> None:
    """Report current (not historical peak) RSS for memory ablations."""
    if os.environ.get("REVAC_LABRADOR_MEMORY_PROFILE", "0") == "0":
        return
    try:
        with open("/proc/self/statm", "r", encoding="ascii") as statm:
            resident_pages = int(statm.read().split()[1])
        rss_mib = resident_pages * os.sysconf("SC_PAGE_SIZE") / (1024 * 1024)
        print(f"[Labrador worker 内存] {label}：当前 RSS={rss_mib:.2f} MiB",
              file=sys.stderr, flush=True)
    except (OSError, ValueError, IndexError):
        # Memory diagnostics must never alter proof availability.
        pass


@contextlib.contextmanager
def polxvec_parallel_scope(phase: str):
    """Select the native vector kernel policy for one worker request.

    Verification operates on immutable proof data and benefits from the long
    challenge-times-vector kernel.  Proving keeps it off because its trial
    loops already expose other teams.  Stress tests cover the verifier path
    with the non-reentrant Gaussian sampler kept serial.
    """
    variable = f"REVAC_LABRADOR_POLXVEC_PARALLEL_{phase.upper()}"
    default = "1" if phase == "verify" else "0"
    enabled = os.environ.get(variable, default)
    previous = os.environ.get("REVAC_LABRADOR_POLXVEC_PARALLEL")
    os.environ["REVAC_LABRADOR_POLXVEC_PARALLEL"] = enabled
    try:
        yield
    finally:
        if previous is None:
            os.environ.pop("REVAC_LABRADOR_POLXVEC_PARALLEL", None)
        else:
            os.environ["REVAC_LABRADOR_POLXVEC_PARALLEL"] = previous


@contextlib.contextmanager
def verify_aggregation_parallel_scope(phase: str):
    """Enable disjoint-output aggregation kernels only for verification.

    Rotation blocks and CRT zeroing write disjoint outputs.  Keep the prover
    serial, and enable the verifier after stress testing it independently from
    the non-reentrant Gaussian sampler.  Phase-specific variables retain an
    immediate rollback and diagnostic knob.
    """
    variables = (
        "REVAC_LABRADOR_ROTATION_PARALLEL",
        "REVAC_LABRADOR_ZERO_PARALLEL",
    )
    previous = {name: os.environ.get(name) for name in variables}
    default = "1" if phase == "verify" else "0"
    for name in variables:
        phase_name = f"{name}_{phase.upper()}"
        os.environ[name] = os.environ.get(phase_name, default)
    try:
        yield
    finally:
        for name in variables:
            old = previous[name]
            if old is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = old


@contextlib.contextmanager
def native_thread_scope(phase: str):
    """Apply a phase-specific native OpenMP budget without changing proofs.

    Proving already runs concurrently with the anonymous-signature prover and
    keeps LaBRADOR's established twelve-thread process default.  Verification
    is the critical path at the GCS after the signature branch became faster.
    A wider 24/16 team improved isolated Pack timings but regressed the real
    complete-show benchmark through contention and longer tails.  Explicit
    phase/global settings always take precedence, which also gives other hosts
    a deterministic tuning and rollback knob.
    """
    phase_variable = f"REVAC_LABRADOR_THREADS_{phase.upper()}"
    phase_value = os.environ.get(phase_variable)
    global_value = os.environ.get("REVAC_LABRADOR_THREADS")
    previous = global_value

    if phase_value is not None:
        selected = phase_value
    elif global_value is not None:
        selected = global_value
    elif phase == "verify":
        selected = str(min(8, max(1, os.cpu_count() or 1)))
    else:
        selected = None

    if selected is not None:
        os.environ["REVAC_LABRADOR_THREADS"] = selected
    try:
        yield
    finally:
        if previous is None:
            os.environ.pop("REVAC_LABRADOR_THREADS", None)
        else:
            os.environ["REVAC_LABRADOR_THREADS"] = previous


def process(mode: bytes, dynamic: bytes, entry: dict[str, Any]) -> tuple[bytes, bool]:
    decode_started = time.perf_counter()
    data, proof_in = decode_dynamic(dynamic, entry["public"], mode)
    if os.environ.get("REVAC_LABRADOR_PROFILE", "0") != "0":
        print(
            "[Labrador worker 计时] 动态请求解码："
            f"{time.perf_counter() - decode_started:.6f} 秒",
            file=sys.stderr, flush=True)
    # The imported implementation is verbose for command-line diagnostics.
    # The persistent service deliberately suppresses that output; C prints one
    # concise bridge timing line per request instead.
    diagnostics = io.StringIO()
    with contextlib.redirect_stdout(diagnostics):
        if mode == b"P":
            proof_out: list[bytes] = []
            generated_params: list[LiveLabrador38Params] = []
            params_cache: OrderedDict[
                tuple[int, int, int, bool, str], LiveLabrador38Params
            ] = entry["params"]
            cache_key = parameter_cache_key(data)
            cached_params = params_cache.get(cache_key)
            if cached_params is not None:
                params_cache.move_to_end(cache_key)
            build_started = time.perf_counter()
            with (native_thread_scope("prove"),
                  polxvec_parallel_scope("prove"),
                  verify_aggregation_parallel_scope("prove")):
                ok = build_labrador(data, data["lift"], self_verify=False,
                                    right_layout="grouped", prover_zk=True,
                                    proof_bytes_out=proof_out,
                                    static_template=ensure_template(entry),
                                    cached_params_obj=cached_params,
                                    params_obj_out=generated_params,
                                    reuse_comkey=True,
                                    compact_params=True)
            if os.environ.get("REVAC_LABRADOR_PROFILE", "0") != "0":
                captured = diagnostics.getvalue()
                if captured:
                    print(captured, file=sys.stderr, end="", flush=True)
                print(
                    "[Labrador worker 计时] statement 绑定 + Pack Prove + 序列化："
                    f"{time.perf_counter() - build_started:.6f} 秒",
                    file=sys.stderr, flush=True)
            if generated_params:
                previous = params_cache.get(cache_key)
                if previous is not None:
                    generated_params[0].close()
                    raise ValueError("同一公开参数策略重复生成了 Labrador 参数")
                params_cache[cache_key] = generated_params[0]
                params_cache.move_to_end(cache_key)
                while len(params_cache) > PARAM_CACHE_LIMIT:
                    _, evicted = params_cache.popitem(last=False)
                    evicted.close()
            if not ok or len(proof_out) != 1:
                raise ValueError("Labrador prove failed")
            entry["active_statement_profile"] = cache_key
            compact_public_after_freeze(entry)
            if os.environ.get("REVAC_LABRADOR_PROFILE", "0") != "0":
                compact = labrador38_proof_is_compact(proof_out[0])
                params = (b"" if compact else
                          labrador38_proof_params_bytes(proof_out[0]))
                params_len = len(params)
                envelope_overhead = len(proof_out[0]) - params_len
                packed_total = len(zlib.compress(proof_out[0], level=9))
                packed_params = (0 if compact else
                                 len(zlib.compress(params, level=9)))
                # The envelope contains magic, two u64 lengths and the inner
                # proof.  Report the split so wire-size optimizations are
                # evaluated against measured bytes rather than estimates.
                print(
                    "[Labrador worker 大小] "
                    f"模式={'compact-v3' if compact else 'full-v2'}；"
                    f"总 envelope={len(proof_out[0])} 字节；"
                    f"公共参数={params_len} 字节；"
                    f"proof+封装={envelope_overhead} 字节；"
                    f"deflate 总计={packed_total} 字节；"
                    f"公共参数单独 deflate={packed_params} 字节",
                    file=sys.stderr,
                    flush=True,
                )
            return proof_out[0], cached_params is not None
    raise ValueError("unsupported MemVer worker mode")


def process_verify(mode: bytes, dynamic: bytes, entry: dict[str, Any]) -> bool:
    """Verify the received binary envelope without temporary files."""
    data, proof_bytes = decode_dynamic(dynamic, entry["public"], mode)
    assert proof_bytes is not None
    # Reject an incompatible compact profile before constructing the 1661
    # witness-polynomial slots and 36 constraints.  The full verifier repeats
    # this check, so this is an early-rejection optimization rather than a
    # replacement for cryptographic verification.
    if labrador38_proof_is_compact(proof_bytes):
        if labrador38_proof_profile_id(proof_bytes) != \
                memver_labrador_profile_id(data):
            return False
    params_cache: OrderedDict[
        tuple[int, int, int, bool, str], LiveLabrador38Params
    ] = entry["params"]
    cache_key = parameter_cache_key(data)
    cached_params = params_cache.get(cache_key)
    if os.environ.get("REVAC_LABRADOR_VERIFY_FORCE_PARAM_MISS", "0") != "0":
        forced = params_cache.pop(cache_key, None)
        if forced is not None:
            forced.close()
        cached_params = None
    if cached_params is None:
        if not labrador38_proof_is_compact(proof_bytes):
            cached_params = LiveLabrador38Params.from_bytes(
                labrador38_proof_params_bytes(proof_bytes),
                require_zk=True, reset_comkey=False)
            params_cache[cache_key] = cached_params
            while len(params_cache) > PARAM_CACHE_LIMIT:
                _, evicted = params_cache.popitem(last=False)
                evicted.close()
    else:
        params_cache.move_to_end(cache_key)
    generated_params: list[LiveLabrador38Params] = []
    if (os.environ.get("REVAC_LABRADOR_FAST_VERIFY", "1") != "0" and
            cached_params is not None and
            entry.get("active_statement_profile") == cache_key):
        fast_started = time.perf_counter()
        with (native_thread_scope("verify"),
              polxvec_parallel_scope("verify"),
              verify_aggregation_parallel_scope("verify")):
            fast_result = verify_labrador_live_statement(
                data, proof_bytes, ensure_template(entry), cached_params,
                reuse_comkey=True)
        if fast_result is not None:
            if os.environ.get("REVAC_LABRADOR_PROFILE", "0") != "0":
                print(
                    "[Labrador worker 计时] 常驻 statement 快速 verify："
                    f"{time.perf_counter() - fast_started:.6f} 秒",
                    file=sys.stderr, flush=True)
            return fast_result
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            with (native_thread_scope("verify"),
                  polxvec_parallel_scope("verify"),
                  verify_aggregation_parallel_scope("verify")):
                ok = bool(build_labrador(
                    data, None, proof_bytes_in=proof_bytes,
                    right_layout="grouped", prover_zk=True,
                    static_template=ensure_template(entry),
                    cached_params_obj=cached_params,
                    params_obj_out=generated_params,
                    reuse_comkey=True))
    except Exception as exc:
        if os.environ.get("REVAC_LABRADOR_PROFILE", "0") != "0":
            print(f"[Labrador worker verify 异常] {exc!r}",
                  file=sys.stderr, flush=True)
        raise
    if os.environ.get("REVAC_LABRADOR_PROFILE", "0") != "0":
        print(f"[Labrador worker verify 结果] {'通过' if ok else '失败'}",
              file=sys.stderr, flush=True)
    if generated_params:
        if not ok or cached_params is not None:
            generated_params[0].close()
            if cached_params is not None:
                raise ValueError("compact verify 重复生成了参数缓存")
        else:
            params_cache[cache_key] = generated_params[0]
            params_cache.move_to_end(cache_key)
            while len(params_cache) > PARAM_CACHE_LIMIT:
                _, evicted = params_cache.popitem(last=False)
                evicted.close()
    compact_public_after_freeze(entry)
    if ok:
        entry["active_statement_profile"] = cache_key
    return ok


def serve(socket_path: str, cache_limit: int) -> int:
    cache: OrderedDict[bytes, dict[str, Any]] = OrderedDict()
    if os.path.exists(socket_path):
        try:
            probe = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            probe.connect(socket_path)
        except OSError:
            os.unlink(socket_path)
        else:
            probe.close()
            raise RuntimeError(f"worker already owns {socket_path}")
    server = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    server.bind(socket_path)
    os.chmod(socket_path, 0o600)
    server.listen(16)
    idle_trim_seconds = float(os.environ.get(
        "REVAC_LABRADOR_IDLE_TRIM_SECONDS", "0.25"))
    if not 0.01 <= idle_trim_seconds <= 60.0:
        raise ValueError(
            "REVAC_LABRADOR_IDLE_TRIM_SECONDS must be in [0.01, 60]")
    trim_policy = os.environ.get(
        "REVAC_LABRADOR_TRIM_POLICY", "idle").strip().lower()
    if trim_policy not in ("idle", "immediate", "retain", "burst"):
        raise ValueError(
            "REVAC_LABRADOR_TRIM_POLICY must be idle, immediate, retain, "
            "or burst")
    cache_idle_seconds = float(os.environ.get(
        "REVAC_LABRADOR_CACHE_IDLE_SECONDS", "30"))
    if not 1.0 <= cache_idle_seconds <= 3600.0:
        raise ValueError(
            "REVAC_LABRADOR_CACHE_IDLE_SECONDS must be in [1, 3600]")
    heap_dirty = False
    last_activity = time.monotonic()

    def cleanup(_signum: int, _frame: Any) -> None:
        for entry in cache.values():
            close_entry(entry)
        cache.clear()
        reset_global_comkey()
        server.close()
        try:
            os.unlink(socket_path)
        except FileNotFoundError:
            pass
        raise SystemExit(0)

    signal.signal(signal.SIGTERM, cleanup)
    signal.signal(signal.SIGINT, cleanup)
    while True:
        # Prove is normally followed immediately by verify.  Returning large
        # temporary arenas to glibc between those two requests used to put
        # gc.collect()+malloc_trim directly on the end-to-end show latency.
        # Keep the same reclamation policy, but perform it only after the
        # socket has remained idle.  A queued request always wins over trim.
        wait_timeout: float | None = None
        now = time.monotonic()
        if heap_dirty and trim_policy in ("idle", "burst"):
            wait_timeout = max(
                0.0, idle_trim_seconds - (now - last_activity))
        if trim_policy == "burst" and cache:
            cache_wait = max(
                0.0, cache_idle_seconds - (now - last_activity))
            wait_timeout = (cache_wait if wait_timeout is None else
                            min(wait_timeout, cache_wait))
        if wait_timeout is not None:
            readable, _, _ = select.select([server], [], [], wait_timeout)
            if not readable:
                idle_for = time.monotonic() - last_activity
                if (trim_policy == "burst" and cache and
                        idle_for >= cache_idle_seconds):
                    for entry in cache.values():
                        close_entry(entry)
                    cache.clear()
                    reset_global_comkey()
                    heap_dirty = False
                    trim_worker_heap()
                    report_worker_memory("长空闲期缓存释放后")
                elif heap_dirty and idle_for >= idle_trim_seconds:
                    trim_worker_heap()
                    report_worker_memory("空闲期堆回收后")
                    heap_dirty = False
                continue
        conn, _ = server.accept()
        with conn:
            try:
                request = recv_frame(conn)
                if len(request) < 1 + HASH_BYTES + 4:
                    raise ValueError("truncated MemVer worker request")
                command, pp_id = request[:1], request[1:1 + HASH_BYTES]
                dynamic_len = struct.unpack_from("<I", request, 1 + HASH_BYTES)[0]
                off = 1 + HASH_BYTES + 4
                if dynamic_len == 0 or off + dynamic_len > len(request):
                    raise ValueError("invalid MemVer worker dynamic length")
                dynamic = request[off:off + dynamic_len]
                entry = cache.get(pp_id)
                if command in (b"P", b"V") and entry is None:
                    send_frame(conn, b"N")
                    continue
                if command in (b"I", b"J"):
                    public_len_off = off + dynamic_len
                    if public_len_off + 4 > len(request):
                        raise ValueError("missing MemVer public-parameter length")
                    public_len = struct.unpack_from("<I", request, public_len_off)[0]
                    public_start = public_len_off + 4
                    if public_len == 0 or public_start + public_len != len(request):
                        raise ValueError("invalid MemVer public-parameter length")
                    entry = {
                        "public": decode_public(request[public_start:public_start + public_len]),
                        "template": None,
                        "params": OrderedDict(),
                        "active_statement_profile": None,
                    }
                    old_entry = cache.get(pp_id)
                    if old_entry is not None:
                        close_entry(old_entry)
                    cache[pp_id] = entry
                    cache.move_to_end(pp_id)
                    while len(cache) > cache_limit:
                        _, evicted = cache.popitem(last=False)
                        close_entry(evicted)
                    command = b"P" if command == b"I" else b"V"
                assert entry is not None
                cache.move_to_end(pp_id)
                if command == b"P":
                    proof, params_hit = process(command, dynamic, entry)
                    send_frame(conn, b"S" + bytes((1 if params_hit else 0,)) + proof)
                    release_parameter_lease(entry)
                    if trim_policy == "immediate":
                        trim_worker_heap()
                        report_worker_memory("prove 请求完成后")
                    else:
                        heap_dirty = True
                elif command == b"V":
                    send_frame(conn, b"S\x00" + (b"\x01" if process_verify(command, dynamic, entry) else b"\x00"))
                    release_parameter_lease(entry)
                    if trim_policy == "immediate":
                        trim_worker_heap()
                        report_worker_memory("verify 请求完成后")
                    else:
                        heap_dirty = True
                else:
                    raise ValueError("unknown MemVer worker command")
            except Exception as exc:  # Return an explicit failure to C.
                send_frame(conn, b"F" + str(exc).encode("utf-8", "replace")[:1024])
        last_activity = time.monotonic()


def main() -> int:
    # A few native helper loops intentionally use OpenMP's default team.  On
    # 64/128-core hosts an unrestricted team made the small refresh-witness
    # composition ~70x slower through launch/NUMA overhead.  Keep those small
    # helper loops and Pack Prove at twelve.  Pack Verify receives its measured
    # eight-thread budget in native_thread_scope().
    os.environ.setdefault(
        "OMP_NUM_THREADS", str(min(12, max(1, os.cpu_count() or 1))))
    parser = argparse.ArgumentParser()
    parser.add_argument("--socket", required=True)
    parser.add_argument("--cache-limit", type=int,
                        default=int(os.environ.get("REVAC_LABRADOR_WORKER_CACHE", "4")))
    args = parser.parse_args()
    if args.cache_limit < 1 or args.cache_limit > 16:
        raise ValueError("--cache-limit must be between 1 and 16")
    # Experimental only: changing glibc allocation layout exposed an upstream
    # LaBRADOR-tail lifetime/uninitialized-state bug in randomized trials.
    # Keep the hook for diagnosis but never enable it in the safe default.
    if os.environ.get("REVAC_LABRADOR_MALLOC_TUNING", "0") != "0":
        from lazer import lib
        result = lib.labrador38_py_malloc_configure(
            int(os.environ.get("REVAC_LABRADOR_MMAP_THRESHOLD", "32768")),
            int(os.environ.get("REVAC_LABRADOR_TRIM_THRESHOLD", "65536")),
            int(os.environ.get("REVAC_LABRADOR_ARENA_MAX", "1")))
        if result != 0:
            raise RuntimeError("无法配置 Labrador worker 的 glibc allocator")
    return serve(args.socket, args.cache_limit)


if __name__ == "__main__":
    raise SystemExit(main())
