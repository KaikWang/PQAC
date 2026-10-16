"""Public production defaults for the complete q88 compact CRT profile.

Set REVAC_Q88_OPTIMIZED_PROFILE=0 in a clean environment to select the
previous compiler layout. Prover and verifier must use the same layout.
All layout changes are bound into the worker static_profile_id.
"""
import os

DEFAULTS = {
    "REVAC_Q88_TIGHT_WITNESS_TOP": "1",
    "REVAC_Q96_AGGREGATE_TIGHT_NORM": "1",
    "REVAC_Q88_NORM_GROUP_SIZE": "16",
    "REVAC_Q88_APPROX_ORDER": "limb",
    "REVAC_Q88_SCALAR_HANDLE_NORM": "1",
    "REVAC_Q88_TIGHT_REFRESH_NORM": "1",
    "REVAC_Q88_AGGREGATE_RADIX": "3779",
    "REVAC_Q88_PRESERVE_AUX_REQUIREMENT": "1",
    "REVAC_Q88_AUX_TOP_UNBALANCED": "1",
    "REVAC_Q88_SCALED_REFRESH_Z": "1",
    "REVAC_Q88_BATCH_SPLIT": "1",
    "REVAC_Q88_FUSED_LING": "1",
    "REVAC_Q88_FUSED_LING_THREADS": "4",
    "REVAC_Q88_SPROD_THREADS": "6",
    "REVAC_LABRADOR_LINEAR_THREADS": "4",
    "REVAC_LABRADOR_SCRATCH_POOL_MIB": "1024",
    "REVAC_LABRADOR_LDD_ADDCHECK_THREADS": "1",
}

def configure():
    if (os.environ.get("REVAC_Q88_OPTIMIZED_PROFILE", "1") == "0" or
            os.environ.get("REVAC_Q96_FULL_RING", "1") == "0"):
        return
    for key, value in DEFAULTS.items():
        os.environ.setdefault(key, value)
