#!/usr/bin/env python3
"""<= 60s Deterministic Verification Harness for Epigenetic State Dispatch v1.0.0."""
import subprocess, sys, time
from pathlib import Path

t0 = time.perf_counter()
here = Path(__file__).resolve().parent
subprocess.run(
    [sys.executable, str(here / "methylated_interior_search_benchmark_v1_0_0.py"),
     "--seed", "42", "--suite", "exp3", "--trials", "12", "--out", str(here / "verify_exp3_smoke.json")],
    check=True
)
elapsed = time.perf_counter() - t0
assert elapsed <= 60.0, f"Exceeded 60s ceiling: {elapsed:.2f}s"
print(f"[PASS] Epigenetic State Dispatch v1.0.0 verified in {elapsed:.2f} s (<= 60 s contract).")
