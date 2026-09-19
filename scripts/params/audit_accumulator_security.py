#!/usr/bin/env python3
"""Fail-closed audit for a theorem-backed revocation-accumulator profile.

The formulas in this file reproduce the *standard M-SIS* branch of
Theorems 4.2 and 4.3 of Kemmoe--Lysyanskaya (PKC 2026).  They are useful for
rejecting impossible parameter proposals and for locating candidates worth a
full estimator run.  Passing this script is deliberately not called a security
certification: the project's balanced gadget and truncated sampler still need
an explicit theorem-to-implementation argument, and the two M-LWE instances
must be evaluated separately.

Notation follows the paper, generalized to two explicitly separated gadgets:
the balanced evaluation gadget controls ``m_prime`` and witness updates, while
the Jeudy--Sanders truncated gadget controls the left trapdoor width.  This
separation matches the implementation, but remains a project lemma rather than
a verbatim instantiation of PKC'26 Theorem 4.3; consequently the audit stays
fail-closed.  ``ring_degree`` is PARAM_N and ``module_rows`` is PARAM_D.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


def ceil_log(value: int, base: int) -> int:
    power = 1
    digits = 0
    while power < value:
        power *= base
        digits += 1
    return digits


def is_prime_64(n: int) -> bool:
    """Deterministic Miller--Rabin for unsigned 64-bit integers."""
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for prime in small:
        if n % prime == 0:
            return n == prime
    d = n - 1
    shift = 0
    while d % 2 == 0:
        d //= 2
        shift += 1
    # Jim Sinclair's deterministic base set for n < 2^64.
    for base in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if base % n == 0:
            continue
        x = pow(base, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(shift - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True


def ntt_prime_below(bits: int, ring_degree: int) -> int:
    """Largest prime q < 2^bits with q = 1 mod 2*ring_degree."""
    modulus = 2 * ring_degree
    candidate = ((1 << bits) - 1) // modulus * modulus + 1
    while candidate.bit_length() == bits:
        if is_prime_64(candidate):
            return candidate
        candidate -= modulus
    raise ValueError(f"no {bits}-bit NTT prime found")


def balanced_decompose(value: int, base: int, digits: int) -> tuple[int, ...]:
    """Canonical centered decomposition used by the proposed odd-base G^-1."""
    if base < 3 or base % 2 == 0:
        raise ValueError("balanced base must be odd and at least three")
    bound = (base - 1) // 2
    out: list[int] = []
    remaining = value
    for _ in range(digits):
        digit = remaining % base
        if digit > bound:
            digit -= base
        out.append(digit)
        remaining = (remaining - digit) // base
    if remaining != 0:
        raise ValueError("insufficient balanced-gadget dimension")
    return tuple(out)


def balanced_kat(q: int, base: int, digits: int) -> bool:
    """Boundary and deterministic pseudo-random known-answer checks."""
    half = (q - 1) // 2
    probes = [-half, -half + 1, -1, 0, 1, half - 1, half]
    state = 0x9E3779B97F4A7C15
    for _ in range(1024):
        state = (state * 6364136223846793005 + 1442695040888963407) \
            & ((1 << 128) - 1)
        probes.append(int(state % q) - half)
    bound = (base - 1) // 2
    for value in probes:
        decomposition = balanced_decompose(value, base, digits)
        if any(abs(digit) > bound for digit in decomposition):
            return False
        if sum(digit * base**i for i, digit in enumerate(decomposition)) != value:
            return False
    return True


def rhf_from_bkz(blocksize: int) -> float:
    b = float(blocksize)
    return (b / (2 * math.pi * math.e) * (math.pi * b) ** (1 / b)) ** (
        1 / (2 * (b - 1))
    )


def find_bkz_blocksize(delta: float) -> int:
    if delta < 1.00057456569590:
        return 4935
    blocksize = 36
    while rhf_from_bkz(blocksize + 100) > delta:
        blocksize += 100
    while rhf_from_bkz(blocksize + 10) > delta:
        blocksize += 10
    while rhf_from_bkz(blocksize + 2) > delta:
        blocksize += 2
    while rhf_from_bkz(blocksize + 1) > delta:
        blocksize += 1
    return blocksize + 1


@dataclass(frozen=True)
class Audit:
    q: int
    q_bits: int
    q_factors: tuple[int, ...]
    ring_degree: int
    module_rows: int
    deletion_bound: int
    eval_base: int
    eval_digit_bound: int
    trunc_base: int
    trapdoor_spectral_slack: float
    eval_digits: int
    trunc_digits: int
    trapdoor_cols: int
    beta0_log2: float
    beta1_log2: float
    beta1_margin_log2: float
    beta1_below_q: bool
    bkz_blocksize: int
    classical_core_svp_bits: float
    quantum_core_svp_bits: float
    rank_failure_bits: float
    accumulator_kib: float
    full_witness_mib: float
    ba_mib: float
    relation_polynomials: int
    balanced_decomposition_kat: bool
    theorem_screen_pass: bool
    release_ready: bool
    blockers: tuple[str, ...]


def audit(
    q: int,
    *,
    ring_degree: int,
    module_rows: int,
    deletion_bound: int,
    handle_bits: int = 32,
    eval_base: int = 3,
    eval_digit_bound: int = 1,
    trunc_base: int = 14,
    trunc_l: int = 2,
    statistical_bits: int = 128,
    spectral_slack: float = 6.0,
    target_quantum_bits: int = 128,
    q_factors: tuple[int, ...] | None = None,
) -> Audit:
    d = ring_degree
    n = module_rows
    ell = handle_bits
    eval_digits = ceil_log(q, eval_base)
    trunc_digits = ceil_log(q, trunc_base)
    if trunc_digits <= trunc_l:
        raise ValueError("truncation requires k > l")
    m_prime = n * eval_digits
    m = n * (2 + trunc_digits - trunc_l)

    # Lemma 4.1 / Theorem 4.2, with the tighter balanced-digit bound.  The
    # substitution eval_digit_bound=1 is conditional on proving the balanced
    # EvalFx analogue of Theorem 2.2, Equation (2.2).
    t = (
        deletion_bound
        * eval_digit_bound
        * math.sqrt(m_prime * d)
        * math.sqrt(
            1
            + (
                eval_digit_bound * d * m_prime * math.sqrt(ell - 1)
            )
            ** 2
        )
    )
    theta1 = math.sqrt(2 * d / 3) * (
        math.sqrt(m) + math.sqrt(m_prime) + math.log(math.sqrt(d))
    )
    theta2 = math.sqrt(2 * d / 3) * (
        math.sqrt(m) + math.sqrt(ell * m_prime) + math.log(math.sqrt(d))
    )
    # Strict lower bound on varsigma is approached from above for screening.
    varsigma = 2 * d * q ** (n / m + 2 / (d * m))
    log_sqrt_nd = math.log(math.sqrt(n * d))

    # RingTruncatedSampler parameters (ePrint 2024/1952, Theorem 4.2),
    # specialized to the accumulator's constant unit tag.  Taking the maximum
    # component width is conservative when screening the norm terms inherited
    # from the spherical SampleLeft presentation.
    eta = math.sqrt(
        (math.log(2 * d * n * trunc_digits)
         + statistical_bits * math.log(2)) / math.pi
    )
    norm_ri = 0.75 * (
        math.sqrt(n * d)
        + math.sqrt((trunc_digits - trunc_l) * n * d)
        + spectral_slack
    )
    trunc_factor = eta * (trunc_base + 1 / trunc_base)
    trunc_s1 = trunc_factor * math.sqrt(4 + 3 * norm_ri * norm_ri)
    trunc_s3 = trunc_factor * math.sqrt(3) * norm_ri
    trunc_s4 = trunc_factor * math.sqrt(3)
    trunc_s11 = trunc_factor * math.sqrt(
        4 * (trunc_base ** (2 * trunc_l) - 1)
        / (trunc_base * trunc_base - 1)
        + 3 * norm_ri * norm_ri
    )
    truncated_width_bound = max(trunc_s1, trunc_s3, trunc_s4, trunc_s11)
    s = max(
        t,
        truncated_width_bound,
        trunc_base * math.sqrt(1 + theta1 * theta1) * log_sqrt_nd,
        trunc_base
        * varsigma
        * math.sqrt(d)
        * (math.sqrt(m) + math.sqrt(m_prime) + math.log(math.sqrt(d)))
        * log_sqrt_nd,
    )
    beta0 = s * math.sqrt(d * (m + ell * m_prime)) + t
    beta1 = (
        (s + 1) * math.sqrt(m * d)
        + theta2 * (s * math.sqrt(d * ell * m_prime) + t)
        + beta0 * math.sqrt(1 + theta2 * theta2)
    )
    beta1_below_q = beta1 < q

    flat_rows = n * d
    flat_cols = m * d
    blocksize = 0
    classical = 0.0
    quantum = 0.0
    if beta1_below_q:
        optimal = round(2 * flat_rows * math.log2(q) / math.log2(beta1))
        optimal = max(flat_rows, min(flat_cols, optimal))
        delta = (beta1 * q ** (-flat_rows / optimal)) ** (1 / optimal)
        blocksize = find_bkz_blocksize(delta)
        classical = 0.292 * blocksize
        quantum = 0.257 * blocksize

    factors = q_factors if q_factors is not None else (q,)
    # For distinct prime factors p_i = 1 mod 2d, every inertia degree is one;
    # Lemma 2.4 gives Pr[rank failure] <= d*sum_i p_i^{-(m-n)}.
    rank_term_logs = tuple(-(m - n) * math.log2(p) for p in factors)
    largest_rank_log = max(rank_term_logs)
    rank_log2 = (
        math.log2(d)
        + largest_rank_log
        + math.log2(sum(2 ** (x - largest_rank_log) for x in rank_term_logs))
    )
    rank_failure_bits = -rank_log2
    q_bits = q.bit_length()
    vector_blocks = 2 + trunc_digits - trunc_l + ell * eval_digits
    accumulator_kib = n * d * q_bits / 8 / 1024
    full_witness_mib = vector_blocks * n * d * q_bits / 8 / 1024 / 1024
    ba_mib = ell * eval_digits * n * n * d * q_bits / 8 / 1024 / 1024
    # Degree-256-equivalent witness polynomials in the current grouped-right
    # compiler.  105 non-witness/static slots are calibrated from the active
    # build; every accumulator polynomial of larger degree occupies d/256
    # backend polynomials.  This is a relation-size projection, not a packed
    # proof-size claim.
    if d % 256:
        raise ValueError("relation projection requires ring degree divisible by 256")
    relation_polynomials = (
        105 + vector_blocks * n * (d // 256)
    )
    decomposition_kat = balanced_kat(q, eval_base, eval_digits)

    blockers: list[str] = []
    if math.prod(factors) != q:
        blockers.append("declared CRT factors do not multiply to q")
    if len(set(factors)) != len(factors):
        blockers.append("q contains repeated prime factors; Lemma 2.4 needs distinct primes")
    if any(p >= (1 << 64) or not is_prime_64(p) for p in factors):
        blockers.append("at least one declared CRT factor is not a verified 64-bit prime")
    if any(p % (2 * d) != 1 for p in factors):
        blockers.append("a CRT factor is not 1 mod 2d; inertia/NTT condition unaudited")
    if not beta1_below_q:
        blockers.append("Theorem 4.3 extraction bound beta1 is not below q")
    if quantum < target_quantum_bits:
        blockers.append(
            f"M-SIS Core-SVP screen is {quantum:.1f} < "
            f"{target_quantum_bits} quantum bits"
        )
    if rank_failure_bits < 128:
        blockers.append("Lemma 2.4 rank-failure probability exceeds 2^-128")
    if not decomposition_kat:
        blockers.append("balanced decomposition known-answer test failed")
    blockers.extend(
        (
            "balanced EvalF/EvalFx norm lemma is not yet connected to the PKC'26 reduction",
            "truncated-sampler covariance and statistical loss are not yet substituted into Theorem 4.3",
            "the separately generated M-LWE rough screen and all reduction/query losses remain to be composed",
            "the integrated same-handle LaZer extraction loss remains to be composed",
        )
    )
    theorem_screen_pass = beta1_below_q and quantum >= target_quantum_bits
    # Fail closed: the analytical integration blockers intentionally prevent
    # code generation from labelling a profile release-ready.
    return Audit(
        q=q,
        q_bits=q_bits,
        q_factors=tuple(factors),
        ring_degree=d,
        module_rows=n,
        deletion_bound=deletion_bound,
        eval_base=eval_base,
        eval_digit_bound=eval_digit_bound,
        trunc_base=trunc_base,
        trapdoor_spectral_slack=spectral_slack,
        eval_digits=eval_digits,
        trunc_digits=trunc_digits,
        trapdoor_cols=m,
        beta0_log2=math.log2(beta0),
        beta1_log2=math.log2(beta1),
        beta1_margin_log2=math.log2(q) - math.log2(beta1),
        beta1_below_q=beta1_below_q,
        bkz_blocksize=blocksize,
        classical_core_svp_bits=classical,
        quantum_core_svp_bits=quantum,
        rank_failure_bits=rank_failure_bits,
        accumulator_kib=accumulator_kib,
        full_witness_mib=full_witness_mib,
        ba_mib=ba_mib,
        relation_polynomials=relation_polynomials,
        balanced_decomposition_kat=decomposition_kat,
        theorem_screen_pass=theorem_screen_pass,
        release_ready=False,
        blockers=tuple(blockers),
    )


def print_audit(label: str, item: Audit) -> None:
    status = "PASS (screen only)" if item.theorem_screen_pass else "REJECT"
    print(f"{label}: {status}")
    print(
        f"  q={item.q} ({item.q_bits} bit), N={item.ring_degree}, "
        f"D={item.module_rows}, gamma={item.deletion_bound}"
    )
    if len(item.q_factors) > 1:
        print("  CRT factors=" + " * ".join(str(p) for p in item.q_factors))
    print(
        f"  eval/trunc digits={item.eval_digits}/{item.trunc_digits}, "
        f"log2(beta0)={item.beta0_log2:.3f}, "
        f"log2(beta1)={item.beta1_log2:.3f}"
    )
    print(f"  log2(q/beta1)={item.beta1_margin_log2:.3f}")
    print(
        f"  M-SIS screen: BKZ={item.bkz_blocksize}, "
        f"classical={item.classical_core_svp_bits:.1f}, "
        f"quantum={item.quantum_core_svp_bits:.1f} bits"
    )
    print(
        f"  structure: accumulator={item.accumulator_kib:.2f} KiB, "
        f"full witness={item.full_witness_mib:.2f} MiB, "
        f"Ba={item.ba_mib:.2f} MiB, "
        f"relation={item.relation_polynomials} polys"
    )
    print(
        "  balanced decomposition KAT="
        + ("pass" if item.balanced_decomposition_kat else "FAIL")
    )
    for blocker in item.blockers:
        print(f"  BLOCKER: {blocker}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gamma", type=int, default=500)
    parser.add_argument("--target-quantum-bits", type=int, default=128)
    parser.add_argument("--json", type=Path)
    parser.add_argument(
        "--scan",
        action="store_true",
        help="scan the current, 24-bit, and conservative PKC-shaped candidates",
    )
    args = parser.parse_args()
    if args.gamma < 1:
        parser.error("--gamma must be positive")

    candidates: list[tuple[str, Audit]] = []
    # Current q is composite and intentionally fails the strict prime profile.
    candidates.append(
        (
            "current 19-bit prototype",
            audit(405769, ring_degree=256, module_rows=4,
                  deletion_bound=args.gamma,
                  target_quantum_bits=args.target_quantum_bits),
        )
    )
    q24 = ntt_prime_below(24, 256)
    candidates.append(
        (
            "24-bit drop-in arithmetic candidate",
            audit(q24, ring_degree=256, module_rows=4,
                  deletion_bound=args.gamma,
                  target_quantum_bits=args.target_quantum_bits),
        )
    )
    if args.scan:
        # The PKC'26 concrete profile uses a rank-one module over a degree-4096
        # ring.  At the same flattened dimension this is preferable here to a
        # smaller ring with larger module rank: Ba is smaller and EvalF/EvalFx
        # perform fewer module-matrix products.  It also minimizes the distance
        # from the theorem's published instantiation.
        for bits in (56, 58, 60, 61, 62):
            q = ntt_prime_below(bits, 4096)
            candidates.append(
                (
                    f"{bits}-bit PKC-shaped conservative screen",
                    audit(q, ring_degree=4096, module_rows=1,
                          deletion_bound=args.gamma,
                          target_quantum_bits=args.target_quantum_bits),
                )
            )
        # A two-limb 80-bit modulus keeps every arithmetic lane below 40 bits,
        # while a large balanced base holds the homomorphic gadget at 8 digits.
        p0 = ntt_prime_below(40, 4096)
        p1 = p0 - 8192
        while not is_prime_64(p1):
            p1 -= 8192
        q80 = p0 * p1
        candidates.append(
            (
                "80-bit two-limb balanced-base-1025 candidate",
                audit(
                    q80,
                    ring_degree=4096,
                    module_rows=1,
                    deletion_bound=args.gamma,
                    eval_base=1025,
                    eval_digit_bound=512,
                    trunc_base=13,
                    q_factors=(p0, p1),
                    target_quantum_bits=args.target_quantum_bits,
                ),
            )
        )
        # Preferred robust profile: three <=32-bit NTT lanes.  Unlike a pair
        # of 48-bit factors, each accumulator congruence can be compiled by
        # the existing 32/36/38-bit LaBRADOR bridge.  The product still covers
        # 96 bits and the odd balanced Eval gadget remains at five digits.
        q96_factors = []
        next_factor = ntt_prime_below(32, 4096)
        while len(q96_factors) < 3:
            q96_factors.append(next_factor)
            next_factor -= 8192
            while not is_prime_64(next_factor):
                next_factor -= 8192
        q96 = math.prod(q96_factors)
        candidates.append(
            (
                "96-bit three-limb balanced-base-602249 candidate",
                audit(
                    q96,
                    ring_degree=4096,
                    module_rows=1,
                    deletion_bound=args.gamma,
                    eval_base=602249,
                    eval_digit_bound=301124,
                    trunc_base=13,
                    spectral_slack=32.0,
                    q_factors=tuple(q96_factors),
                    target_quantum_bits=args.target_quantum_bits,
                ),
            )
        )

    for index, (label, item) in enumerate(candidates):
        if index:
            print()
        print_audit(label, item)

    if args.json is not None:
        payload = {
            "schema_version": 1,
            "warning": "screening result, not a security certification",
            "methodology": {
                "accumulator_reference": "IACR ePrint 2025/1099, Theorems 4.2 and 4.3",
                "core_svp_model": "estimator-style root-Hermite screen; classical 0.292b, quantum 0.257b",
                "balanced_gadget_status": "project lemma; complete reduction integration pending",
            },
            "candidates": [
                {"label": label, **asdict(item)} for label, item in candidates
            ],
        }
        args.json.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
