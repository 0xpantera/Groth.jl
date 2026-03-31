#!/usr/bin/env python3

import importlib.metadata as importlib_metadata
import json
import math
import statistics
import sys
import time
from pathlib import Path


ORDER_R = 21888242871839275222246405745257275088548364400416034343698204186575808495617
ACCUM_SIZES = (32, 128)
TARGET_SECONDS = 0.2
ROUNDS = 7
MAX_LOOPS = 10000


def _patch_py_ecc_version():
    real_version = importlib_metadata.version

    def patched_version(name: str) -> str:
        if name == "py_ecc":
            return "local-checkout"
        return real_version(name)

    importlib_metadata.version = patched_version


def load_py_ecc():
    workspace_root = Path(__file__).resolve().parents[2]
    py_ecc_root = workspace_root / "py_ecc"
    if not py_ecc_root.is_dir():
        raise FileNotFoundError(f"py_ecc checkout not found at {py_ecc_root}")

    _patch_py_ecc_version()
    sys.path.insert(0, str(py_ecc_root))
    from py_ecc.optimized_bn128 import G1, G2, add, multiply, normalize, pairing

    return {
        "G1": G1,
        "G2": G2,
        "add": add,
        "multiply": multiply,
        "normalize": normalize,
        "pairing": pairing,
        "py_ecc_root": str(py_ecc_root),
    }


def det_scalar(tag: int) -> int:
    return ((tag * tag * tag + 19 * tag + 7) % (ORDER_R - 1)) + 1


def make_summary(samples):
    min_seconds = min(samples)
    median_seconds = statistics.median(samples)
    return {
        "min_pretty": f"{min_seconds:.6f} s",
        "median_pretty": f"{median_seconds:.6f} s",
        "min_seconds": min_seconds,
        "median_seconds": median_seconds,
        "memory_bytes": None,
    }


def benchmark(func):
    func()
    t0 = time.perf_counter()
    func()
    first = time.perf_counter() - t0
    loops = max(1, min(MAX_LOOPS, math.ceil(TARGET_SECONDS / max(first, 1e-6))))
    samples = []
    for _ in range(ROUNDS):
        start = time.perf_counter()
        for _ in range(loops):
            func()
        elapsed = time.perf_counter() - start
        samples.append(elapsed / loops)
    return make_summary(samples)


def serialize_g1(normalize, point):
    x, y = normalize(point)
    return [str(x), str(y)]


def serialize_g2(normalize, point):
    x, y = normalize(point)
    return [[str(x.coeffs[0]), str(x.coeffs[1])], [str(y.coeffs[0]), str(y.coeffs[1])]]


def serialize_fq12(element):
    return [str(c) for c in element.coeffs]


def build_accum_inputs(multiply, generator, base_offset: int, scalar_offset: int, size: int):
    bases = [multiply(generator, det_scalar(base_offset + i)) for i in range(1, size + 1)]
    scalars = [det_scalar(scalar_offset + i) for i in range(1, size + 1)]
    return bases, scalars


def naive_accum(add, multiply, bases, scalars):
    acc = None
    for base, scalar in zip(bases, scalars):
        term = multiply(base, scalar)
        acc = term if acc is None else add(acc, term)
    return acc


def is_one_fq12(element) -> bool:
    return element.coeffs[0] == 1 and all(c == 0 for c in element.coeffs[1:])


def main():
    env = load_py_ecc()
    G1 = env["G1"]
    G2 = env["G2"]
    add = env["add"]
    multiply = env["multiply"]
    normalize = env["normalize"]
    pairing = env["pairing"]

    g1_base = multiply(G1, det_scalar(101))
    g2_base = multiply(G2, det_scalar(102))
    g1_scalar = det_scalar(201)
    g2_scalar = det_scalar(202)
    pair_p = multiply(G1, det_scalar(301))
    pair_q = multiply(G2, det_scalar(302))

    results = {
        "scalar_mul": {
            "g1": benchmark(lambda: multiply(g1_base, g1_scalar)),
            "g2": benchmark(lambda: multiply(g2_base, g2_scalar)),
        },
        "naive_accum_g1": {},
        "naive_accum_g2": {},
        "pairing": {
            "single": benchmark(lambda: pairing(pair_q, pair_p)),
        },
    }

    semantic = {
        "scalar_mul": {
            "g1": serialize_g1(normalize, multiply(g1_base, g1_scalar)),
            "g2": serialize_g2(normalize, multiply(g2_base, g2_scalar)),
        },
        "naive_accum_g1": {},
        "naive_accum_g2": {},
        "pairing": {
            "coeffs": serialize_fq12(pairing(pair_q, pair_p)),
        },
        "pairing_checks": {},
    }

    for size in ACCUM_SIZES:
        g1_bases, g1_scalars = build_accum_inputs(multiply, G1, 1000 + size, 2000 + size, size)
        g2_bases, g2_scalars = build_accum_inputs(multiply, G2, 3000 + size, 4000 + size, size)
        g1_acc = naive_accum(add, multiply, g1_bases, g1_scalars)
        g2_acc = naive_accum(add, multiply, g2_bases, g2_scalars)
        results["naive_accum_g1"][str(size)] = benchmark(lambda b=g1_bases, s=g1_scalars: naive_accum(add, multiply, b, s))
        results["naive_accum_g2"][str(size)] = benchmark(lambda b=g2_bases, s=g2_scalars: naive_accum(add, multiply, b, s))
        semantic["naive_accum_g1"][str(size)] = serialize_g1(normalize, g1_acc)
        semantic["naive_accum_g2"][str(size)] = serialize_g2(normalize, g2_acc)

    base_pair = pairing(pair_q, pair_p)
    double_g1 = pairing(pair_q, multiply(pair_p, 2))
    double_g2 = pairing(multiply(pair_q, 2), pair_p)
    semantic["pairing_checks"] = {
        "double_bilinear": double_g1 == double_g2,
        "non_degenerate": not is_one_fq12(base_pair),
    }

    payload = {
        "meta": {
            "python_version": sys.version,
            "py_ecc_root": env["py_ecc_root"],
            "py_ecc_version": "local-checkout",
        },
        "results": results,
        "semantic": semantic,
    }
    json.dump(payload, sys.stdout)


if __name__ == "__main__":
    main()
