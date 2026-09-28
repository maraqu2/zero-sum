Run from this directory with `sage -python sage_data.py --bound <bound> --k <query-budget>`.

For example:

```sh
sage -python sage_data.py --bound 1000 --k 5
sage -python sage_data.py --bound 1000 --k 7
sage -python sage_data.py --bound 1000 --k 11
```

The script examines every odd prime p up to the bound. It asks whether the group C_p of p-th roots of unity in the binary extension field GF(2^ord_p(2)) contains a zero-sum subset of odd size **at most k**. The default is k=3; any odd integer at least 3 is supported, including every odd prime. The field characteristic is two, not p.

Each run creates a timestamped subdirectory of `sage_runs` with:

- `primes.csv`: exact triple counts plus the chosen query budget, `query_exists` (`yes`, `no`, or `unknown`), `Ak` (the count at exactly size k), and `odd_count_le_k` (the total count at odd sizes through k). Blank counts are unknown, never zero. `Ak_status` states whether that count is exact. A known budget answer can coexist with an unknown exact-size count.
- `query_results.json`: the same budget results with the method and field representation for additional enumerations.
- `spectra.json`: full spectra and minimum odd witnesses for selected small primes, together with their query-budget results.
- `metadata.json`: settings, Sage/Python versions, runtime, and counts of yes/no/unknown answers.

Budget k means ell(p)<=k, not A_k(p)>0. For example, at p=7 budget five succeeds through a triple, although A_5(7)=0. At p=23 budget five fails and budget seven succeeds.

The script always computes triples by polynomial gcd. For larger budgets it reuses selected spectra, applies exact elementary identities, or enumerates the smaller of the zero-sum code and its dual. Default enumeration caps are dimension 16 and field degree 20. Unresolved cases are explicitly `unknown`; for example p=337 is unknown at budget five under the default caps. Increasing k does not remove these computational limits.

See [SAGE_WORKFLOW.md](SAGE_WORKFLOW.md) for the methods, options, and validation command, and [SAGE_FLOWCHART.md](SAGE_FLOWCHART.md) for the flow chart.
