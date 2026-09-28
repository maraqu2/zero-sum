"""Exact zero-sum data in characteristic two.

Run: sage -python ldc_zero_sums/sage_data.py --bound 1000 --k 5
The odd query budget k means ell(p) <= k, not necessarily A_k(p) > 0.
Resource-limited results are explicitly unknown; previous runs are preserved.
"""

import argparse
import csv
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter

from sage.all import GF, ZZ, Mod, PolynomialRing, binomial, matrix, prime_range
from sage.env import SAGE_VERSION


def validate_query_k(query_k):
    if query_k < 3 or query_k % 2 == 0:
        raise ValueError("k must be an odd integer at least 3")


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


def field_model(row):
    """Construct the roots and the binary matrix that tests their sums."""
    p, r = row["p"], row["r"]
    F = GF(ZZ(2)**r, name="a")
    zeta = F.multiplicative_generator()**((F.order() - 1) // p)
    roots = [zeta**i for i in range(p)]
    assert zeta != 1 and zeta**p == 1 and len(set(roots)) == p
    # Column i is the coordinate vector of zeta^i over GF(2).
    H = matrix(GF(2), [list(z) for z in roots]).transpose()
    assert H.rank() == r
    description = dict(
        field_modulus_coefficients=[int(a) for a in F.modulus().list()],
        zeta_coordinates=[int(a) for a in zeta],
        coordinate_convention="coefficients in ascending powers of F.gen()")
    return F, roots, H, description


def small_spectrum(p, max_dimension=16, max_field_degree=20):
    row = triple_row(p)
    p, r = row["p"], row["r"]
    dimension = p - r
    if dimension > max_dimension or r > max_field_degree:
        return dict(p=p, r=r, status="not computed: resource cap",
                    dimension=dimension, max_dimension=max_dimension,
                    max_field_degree=max_field_degree)
    started = perf_counter()
    F, roots, H, description = field_model(row)
    kernel = H.right_kernel()
    assert kernel.dimension() == dimension
    A = [0] * (p + 1)
    ell, witness = p + 1, None
    for c in kernel:
        w = int(c.hamming_weight())
        A[w] += 1
        if w % 2 and w < ell:
            ell = w
            witness = [i for i in range(p) if c[i]]
    assert sum(A) == 2**dimension
    assert sum(A[1::2]) == 2**(dimension-1)
    assert A == A[::-1]
    assert A[3] == row["A3"]
    assert A[0] == A[p] == 1 and A[1] == A[2] == 0
    assert witness is not None and sum((roots[i] for i in witness), F(0)) == 0
    return dict(p=p, r=r, dimension=dimension, status="exact exhaustive enumeration",
                **description,
                A=A, odd_spectrum=[s for s in range(3, p+1, 2) if A[s]],
                ell=ell, witness_exponents=witness,
                words_enumerated=sum(A), runtime_seconds=perf_counter()-started)


def counts_through(row, max_size, max_dimension, max_field_degree):
    """Exact A_0,...,A_max_size, enumerating the smaller binary code.

    For dual weight w, the Krawtchouk coefficients K_s(w) satisfy
    (s+1)K_(s+1)=(p-2w)K_s-(p-s+1)K_(s-1). MacWilliams gives
    A_s = 2^(-r) sum_w B_w K_s(w). All arithmetic is integral.
    """
    p, r, dimension = row["p"], row["r"], row["kernel_dimension"]
    limit = min(max_size, p)
    if r > max_field_degree or min(r, dimension) > max_dimension:
        return dict(status="not computed: resource cap",
                    reason=(f"field degree {r} (cap {max_field_degree}); "
                            f"smaller code dimension {min(r, dimension)} "
                            f"(cap {max_dimension})"))
    started = perf_counter()
    _, _, H, description = field_model(row)
    A = [0] * (limit + 1)
    if dimension <= r:
        kernel = H.right_kernel()
        assert kernel.dimension() == dimension
        for c in kernel:
            weight = int(c.hamming_weight())
            if weight <= limit:
                A[weight] += 1
        method = "exhaustive kernel enumeration"
        enumerated_dimension = dimension
    else:
        dual = H.row_space()
        assert dual.dimension() == r
        B = [0] * (p + 1)
        for c in dual:
            B[int(c.hamming_weight())] += 1
        assert sum(B) == 2**r
        numerators = [0] * (limit + 1)
        for weight, multiplicity in enumerate(B):
            if not multiplicity:
                continue
            previous, current = 0, 1
            for s in range(limit + 1):
                numerators[s] += multiplicity * current
                if s < limit:
                    numerator = (p - 2*weight)*current - (p - s + 1)*previous
                    assert numerator % (s + 1) == 0
                    previous, current = current, numerator // (s + 1)
        assert all(n >= 0 and n % (2**r) == 0 for n in numerators)
        A = [n // (2**r) for n in numerators]
        method = "exhaustive dual enumeration and exact MacWilliams transform"
        enumerated_dimension = r
    assert A[0] == 1 and A[1] == A[2] == 0 and A[3] == row["A3"]
    assert all(A[s] <= binomial(p, s) for s in range(limit + 1))
    if limit == p:
        assert A == A[::-1] and sum(A) == 2**dimension
        assert sum(A[1::2]) == 2**(dimension-1)
    return dict(status="exact counts through requested size", method=method,
                A_prefix=A, counted_through=limit, **description,
                enumerated_dimension=enumerated_dimension,
                words_enumerated=2**enumerated_dimension,
                runtime_seconds=perf_counter()-started)


def query_result(row, query_k, max_dimension=16, max_field_degree=20, spectrum=None):
    """Decide the odd size budget separately from the exact-size count.

    None counts mean unknown; zero always means a proved exact zero.
    A positive gcd count can settle existence even when A_k is unknown.
    """
    validate_query_k(query_k)
    p, dimension = row["p"], row["kernel_dimension"]
    result = dict(p=p, query_k=query_k, query_exists="unknown",
                  Ak=None, Ak_status="unknown", odd_count_le_k=None,
                  query_method="not computed: resource cap", query_note="",
                  ell_status=row["ell_status"])
    if spectrum is not None and "A" in spectrum:
        A = spectrum["A"]
        result.update(query_method="exhaustive selected-prime spectrum",
                      ell_status=f"exact: {spectrum['ell']}")
    elif query_k == 3:
        A = [1, 0, 0, row["A3"]]
        result["query_method"] = "exact polynomial gcd over GF(2)"
    elif dimension == 1:
        # Only the empty and full subsets lie in this one-dimensional kernel.
        A = [0] * (min(query_k, p) + 1)
        A[0] = 1
        if query_k >= p:
            A[p] = 1
        result.update(query_method="one-dimensional kernel identity",
                      ell_status=f"exact: {p}")
    else:
        computation = counts_through(row, query_k, max_dimension, max_field_degree)
        result["computation"] = computation
        A = computation.get("A_prefix")
        if A is not None:
            result["query_method"] = computation["method"]
        else:
            result["query_note"] = computation["reason"]
    if A is not None:
        result["Ak"] = A[query_k] if query_k <= p else 0
        result["Ak_status"] = "exact"
        result["odd_count_le_k"] = sum(A[3:min(query_k, p)+1:2])
        result["query_exists"] = "yes" if result["odd_count_le_k"] else "no"
        odd_sizes = [s for s in range(3, min(query_k, p)+1, 2) if A[s]]
        if odd_sizes:
            result["ell_status"] = f"exact: {odd_sizes[0]}"
        elif not result["ell_status"].startswith("exact:"):
            result["ell_status"] = f"lower bound: at least {query_k + 2}"
        return result

    # Elementary certificates still apply when enumeration is unavailable.
    if query_k >= p:
        result.update(query_exists="yes", Ak=int(query_k == p), Ak_status="exact",
                      odd_count_le_k=2**(dimension-1),
                      query_method="whole-group relation and complement identity")
    elif row["A3"]:
        result.update(query_exists="yes", query_method="exact triple within budget")
    if query_k < p and p - query_k <= 2:
        # Complementation reduces this count to nonzero singletons or pairs.
        result.update(Ak=0, Ak_status="exact")
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bound", type=int, default=1000)
    parser.add_argument("--k", type=int, default=3,
                        help="odd query budget >= 3; tests ell(p) <= k (default: 3)")
    parser.add_argument("--small-primes", default="3,5,7,11,17,23")
    parser.add_argument("--max-dimension", type=int, default=16,
                        help="maximum dimension of a code to enumerate (default: 16)")
    parser.add_argument("--max-field-degree", type=int, default=20)
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).resolve().parent / "sage_runs",
                        help="parent directory for new timestamped runs")
    args = parser.parse_args(argv)
    if args.bound < 3:
        parser.error("--bound must be at least 3")
    try:
        validate_query_k(args.k)
        selected = list(dict.fromkeys(int(p) for p in args.small_primes.split(",") if p.strip()))
    except ValueError as exc:
        parser.error(str(exc))
    if any(p < 3 or not ZZ(p).is_prime() for p in selected):
        parser.error("--small-primes must contain only odd primes")
    if args.max_dimension < 0 or args.max_field_degree < 1:
        parser.error("resource caps require max-dimension >= 0 and max-field-degree >= 1")
    started = perf_counter()
    rows = [triple_row(p) for p in prime_range(3, args.bound + 1)]
    samples = [small_spectrum(p, args.max_dimension, args.max_field_degree) for p in selected]
    spectra = {s["p"]: s for s in samples}
    results = [query_result(row, args.k, args.max_dimension, args.max_field_degree,
                            spectra.get(row["p"])) for row in rows]
    for row, result in zip(rows, results):
        row.update({key: value for key, value in result.items()
                    if key not in ("p", "computation")})
    results_by_prime = {result["p"]: result for result in results}
    for sample in samples:
        result = results_by_prime.get(sample["p"])
        if result is None:
            result = query_result(triple_row(sample["p"]), args.k, args.max_dimension,
                                  args.max_field_degree, sample)
        sample["query"] = result
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    out = args.output_dir / stamp
    out.mkdir(parents=True, exist_ok=False)
    with (out / "primes.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (out / "spectra.json").write_text(json.dumps(samples, indent=2) + "\n")
    (out / "query_results.json").write_text(json.dumps(results, indent=2) + "\n")
    summary = {state: sum(row["query_exists"] == state for row in rows)
               for state in ("yes", "no", "unknown")}
    parameters = dict(vars(args), output_dir=str(args.output_dir.resolve()))
    metadata = dict(sage_version=str(SAGE_VERSION), python_version=platform.python_version(),
                    parameters=parameters, timestamp_utc=stamp,
                    query_semantics="an odd zero-sum subset of size at most k",
                    query_summary=summary,
                    sampling="all odd primes <= bound; small spectra separately selected",
                    randomness="no random sampling; field models and generators saved",
                    odd_primes_scanned=len(rows), runtime_seconds=perf_counter()-started)
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    yes = [row["p"] for row in rows if row["A3"]]
    print("Primes with A3 > 0:", yes)
    print("Exact triple fraction:", len(yes), "/", len(rows))
    print(f"Query budget k = {args.k}: " + ", ".join(f"{state}={n}" for state, n in summary.items()))
    if summary["unknown"]:
        print("Unknown cases exceeded resource caps; they are not negative results.")
    for s in samples:
        print("p =", s["p"], ";", s["status"], "; ell =", s.get("ell", "unknown"))
    print("Saved:", out)


if __name__ == "__main__":
    main()
