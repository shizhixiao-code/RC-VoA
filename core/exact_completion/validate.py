"""Validate both directions against hidden-label and literal subset enumeration."""
from pathlib import Path
from itertools import product
from types import SimpleNamespace
from unittest.mock import patch
import json
import numpy as np
import exact_completion as ec

HERE = Path(__file__).resolve().parent


def literal_feasible(labels, cells, ep, bf, fc, K):
    # Excluded endpoints are explicitly represented among labels here.
    b, p, r = ec.B[cells], ec.P[cells], ec.R[cells]
    cancer, benign = (labels >= 0) & (labels < 6), labels == 6
    for bits in product([False, True], repeat=len(labels)):
        if sum(bits) > K: continue
        final = np.where(np.array(bits), p, b)
        if np.any(cancer & (b | r) & ~final): continue
        if final.sum() > fc: continue
        if 10 * (benign & final).sum() > 9 * (benign & b).sum(): continue
        if 10 * (benign & final).sum() > 9 * (benign & r).sum(): continue
        if 100 * (cancer & final).sum() < 90 * cancer.sum(): continue
        if any(100 * ((labels == j) & final).sum() < rate * (labels == j).sum()
               for j, rate in enumerate(ec.RATES)): continue
        return True
    return False


def main():
    rng = np.random.default_rng(20261107)
    scenarios = completions = literals = marginal_checks = 0
    verdict_counts = {-1: 0, 0: 0, 1: 0}
    for rep in range(180):
        N = int(rng.integers(1, 8))
        labels = rng.integers(0, 7, N)
        cells = rng.integers(0, 8, N)
        if rep % 3 == 0: labels[:] = 6; cells[:] = 5
        if rep % 5 == 0: labels[0] = -1
        if rep % 7 == 0: labels[0] = 0; cells[0] = 3
        known = rng.random(N) < 0.6
        known[labels < 0] = True
        known[np.flatnonzero(~known)[3:]] = True
        ix = np.flatnonzero(~known)
        seen = np.zeros((7, 8), np.int64)
        np.add.at(seen, (labels[known & (labels >= 0)], cells[known & (labels >= 0)]), 1)
        unresolved = np.bincount(cells[~known], minlength=8)
        ep = int(((labels < 0) & ec.DOWN[cells]).sum())
        bf = int(ec.B[cells].sum())
        fc = min(int(rng.choice([30, 50, 100])) * N // 100, bf - (N + 9) // 10,
                 int(ec.R[cells].sum()) - (N + 9) // 10)
        K = int(rng.choice([0, N // 5, N // 2, N]))
        margins = []
        for choice in product(range(7), repeat=len(ix)):
            lab = labels.copy(); lab[ix] = choice
            x = np.zeros((7, 8), np.int64)
            valid = lab >= 0
            np.add.at(x, (lab[valid], cells[valid]), 1)
            full = ec.full_margins(x, ep, bf, fc, K)
            if N <= 5:
                assert bool((full >= 0).all()) == literal_feasible(lab, cells, ep, bf, fc, K)
                literals += 1
            margins.append(full)
        margins = np.array(margins)
        truths = (margins >= 0).all(1)
        args = seen, unresolved, ep, bf, fc, K
        status, _, _ = ec.solve_query(*args)
        assert status == ('SAT' if truths.any() else 'UNSAT'), (rep, status)
        for j in range(10):
            status, _, _ = ec.solve_query(*args, margin=j)
            assert status == ('SAT' if (margins[:, j] < 0).any() else 'UNSAT'), (rep, j, status)
            marginal_checks += 1
        result = ec.verify(*args)
        truth = 1 if truths.all() else -1 if not truths.any() else 0
        assert result.verdict == truth and result.unknowns == 0
        verdict_counts[truth] += 1
        scenarios += 1; completions += len(truths)
    # Larger integer ceilings and coupled max branches, fully observed.
    for rep in range(70):
        x = rng.integers(0, 24, (7, 8))
        if rep < 10:
            x[:] = 0; x[6, 5] = [9, 10, 11, 19, 20, 21, 99, 100, 101, 199][rep]
            x[0, 3] = rep
        N = int(x.sum())
        bf, rf = int((x * ec.B).sum()), int((x * ec.R).sum())
        fc = min(N // 2, bf - (N + 9) // 10, rf - (N + 9) // 10)
        K = int(rng.choice([N // 5, N // 2, N]))
        args = x, np.zeros(8, np.int64), 0, bf, fc, K
        actual = ec.full_margins(x, 0, bf, fc, K)
        assert ec.verify(*args).verdict == (1 if (actual >= 0).all() else -1)
        for j in range(10):
            status, _, _ = ec.solve_query(*args, margin=j)
            assert status == ('SAT' if actual[j] < 0 else 'UNSAT')
            marginal_checks += 1
    # Fail-closed solver statuses, including an apparently present incumbent.
    for raw_status in [1, 3, 4]:
        with patch.object(ec, 'milp', return_value=SimpleNamespace(status=raw_status, x=np.zeros(56))):
            result = ec.verify(np.zeros((7, 8), int), np.ones(8, int), 0, 4, 3, 8)
            assert result.verdict == 0 and result.unknowns == 1
    # Also exercise timeout during adverse search after a known SAT pass.
    with patch.object(ec, 'solve_query', side_effect=[('SAT', None, 0), ('UNKNOWN', None, 1)]):
        result = ec.verify(np.zeros((7, 8), int), np.ones(8, int), 0, 4, 3, 8)
        assert result.verdict == 0 and result.unknowns == 1
    report = dict(status='PASS', partial_scenarios=scenarios, completions=completions,
                  literal_allocation_completions=literals, individual_margin_checks=marginal_checks,
                  larger_full_tables=70, verdict_counts=verdict_counts, timeout_and_error_tests=4)
    (HERE / 'VALIDATION.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__': main()
