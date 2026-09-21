import argparse
import csv
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter

from sage.all import GF, ZZ, Mod, PolynomialRing, binomial, matrix, prime_range
from sage.env import SAGE_VERSION


def triple_row(p):
    p = ZZ(p)
    if not p.is_prime() or p == 2:
        raise ValueError("p must be an odd prime")
    r = int(Mod(2, p).multiplicative_order())
    R = PolynomialRing(GF(2), "x")
    x = R.gen()
    g = (x**p + 1).gcd((x + 1)**p + 1)
    D = int(g.degree())
    assert (int(p) * D) % 6 == 0
    A3 = int(p) * D // 6
    return dict(p=int(p), r=r, kernel_dimension=int(p)-r,
                gcd_degree=D, A3=A3,
                triple_probability=str(ZZ(A3) / binomial(p, 3)),
                r_even=(r % 2 == 0),
                lemma8_hypothesis=bool(ZZ(2)**(3*r) < p**4),
                lemma23_necessary_condition=bool(3*r*r <= 4*p),
                ell_status="exact: 3" if A3 else "lower bound: at least 5",
                method="exact polynomial gcd over GF(2)")


def small_spectrum(p, max_dimension=16, max_field_degree=20):
    row = triple_row(p)
    p, r = row["p"], row["r"]
    k = p - r
    if k > max_dimension or r > max_field_degree:
        return dict(p=p, r=r, status="not computed: resource cap",
                    dimension=k, max_dimension=max_dimension,
                    max_field_degree=max_field_degree)
    started = perf_counter()
    F = GF(ZZ(2)**r, name="a")
    zeta = F.multiplicative_generator()**((F.order() - 1) // p)
    roots = [zeta**i for i in range(p)]
    assert zeta != 1 and zeta**p == 1 and len(set(roots)) == p
    # Column i is the coordinate vector of zeta^i over GF(2).
    H = matrix(GF(2), [list(z) for z in roots]).transpose()
    assert H.rank() == r
    kernel = H.right_kernel()
    assert kernel.dimension() == k
    A = [0] * (p + 1)
    ell, witness = p + 1, None
    for c in kernel:
        w = int(c.hamming_weight())
        A[w] += 1
        if w % 2 and w < ell:
            ell = w
            witness = [i for i in range(p) if c[i]]
    assert sum(A) == 2**k
    assert sum(A[1::2]) == 2**(k-1)
    assert A == A[::-1]
    assert A[3] == row["A3"]
    assert A[0] == A[p] == 1 and A[1] == A[2] == 0
    assert witness is not None and sum((roots[i] for i in witness), F(0)) == 0
    return dict(p=p, r=r, dimension=k, status="exact exhaustive enumeration",
                field_modulus_coefficients=[int(a) for a in F.modulus().list()],
                zeta_coordinates=[int(a) for a in zeta],
                coordinate_convention="coefficients in ascending powers of F.gen()",
                A=A, odd_spectrum=[s for s in range(3, p+1, 2) if A[s]],
                ell=ell, witness_exponents=witness,
                words_enumerated=sum(A), runtime_seconds=perf_counter()-started)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bound", type=int, default=1000)
    parser.add_argument("--small-primes", default="3,5,7,11,17,23")
    parser.add_argument("--max-dimension", type=int, default=16)
    parser.add_argument("--max-field-degree", type=int, default=20)
    args = parser.parse_args()
    if args.bound < 3:
        parser.error("--bound must be at least 3")
    started = perf_counter()
    rows = [triple_row(p) for p in prime_range(3, args.bound + 1)]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    out = Path(__file__).resolve().parent / "sage_runs" / stamp
    out.mkdir(parents=True, exist_ok=False)
    with (out / "primes.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    samples = [small_spectrum(int(p), args.max_dimension, args.max_field_degree)
               for p in args.small_primes.split(",") if p.strip()]
    (out / "spectra.json").write_text(json.dumps(samples, indent=2) + "\n")
    metadata = dict(sage_version=str(SAGE_VERSION), python_version=platform.python_version(),
                    parameters=vars(args), timestamp_utc=stamp,
                    sampling="all odd primes <= bound; small spectra separately selected",
                    randomness="no random sampling; field models and generators saved",
                    odd_primes_scanned=len(rows), runtime_seconds=perf_counter()-started)
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    yes = [row["p"] for row in rows if row["A3"]]
    print("Primes with A3 > 0:", yes)
    print("Exact finite-sample fraction:", len(yes), "/", len(rows))
    for s in samples:
        print("p =", s["p"], ";", s["status"], "; ell =", s.get("ell", "unknown"))
    print("Saved:", out)


if __name__ == "__main__":
    main()
