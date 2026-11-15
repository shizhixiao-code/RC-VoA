"""Offline earliest-prefix audit and separately timed, sequential online replay."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
from time import perf_counter
import numpy as np
import pandas as pd
import scipy
import exact_completion as ec
from rc_interval import certificate

HERE = Path(__file__).resolve().parent


def load(key):
    with np.load(HERE / 'inputs' / f'{key}.npz', allow_pickle=False) as z:
        return {k: z[k].copy() for k in z.files}


def exact_args(data, seen):
    N, E, K, ep, bf, fc = map(int, data['public'])
    return seen, data['cell_totals'] - seen.sum(0), ep, bf, fc, K


def offline(meta, order_ids, dest):
    data = load(meta['table_id'])
    N, E, K, ep, bf, fc = map(int, data['public'])
    # Independent full-record solve; not supplied to the prefix search.
    full = ec.verify(*exact_args(data, data['group_counts'].sum(0)))
    assert full.verdict == int(data['truth']) and not full.unknowns
    rows, traces, boundaries = [], [], []
    cache = {}
    start = perf_counter()
    for oi in order_ids:
        order = data['orders'][oi]
        seen = np.concatenate([np.zeros((1, 7, 8), dtype=np.int64), data['group_counts'][order].cumsum(0)])
        costs = seen.sum((1, 2))
        calls = unknowns = evaluations = 0
        seconds = 0.0
        def evaluate(prefix):
            nonlocal calls, unknowns, evaluations, seconds
            key = seen[prefix].tobytes()
            hit = key in cache
            if hit: result = cache[key]
            else:
                result = ec.verify(*exact_args(data, seen[prefix]))
                cache[key] = result
                calls += result.calls; unknowns += result.unknowns; seconds += result.seconds
            evaluations += 1
            traces.append(dict(order=int(oi), prefix=int(prefix), verdict=result.verdict,
                               pass_status=result.possible_pass, adverse_margin=result.adverse_margin,
                               calls=0 if hit else result.calls, unknowns=result.unknowns, cache_hit=hit,
                               seconds=0 if hit else result.seconds))
            # UNKNOWN destroys the equivalence of unresolved and two witnesses.
            # Never silently binary-search through a solver timeout.
            if result.unknowns:
                raise RuntimeError(f'Offline exact bracket unproved: {meta["table_id"]}, order={oi}, prefix={prefix}')
            assert not result.verdict or result.verdict == int(data['truth'])
            if result.verdict == 0:
                assert result.pass_witness is not None and result.adverse_witness is not None
            return result
        if meta['screen']:
            # Shared decision-only pre-screen, verified against a full MILP above.
            stop, result, previous = 0, full, None
        else:
            zero = evaluate(0)
            if zero.verdict:
                stop, result, previous = 0, zero, None
            else:
                lo, hi = 0, len(order)
                assert evaluate(hi).verdict != 0
                while hi - lo > 1:
                    mid = (lo + hi) // 2
                    if evaluate(mid).verdict: hi = mid
                    else: lo = mid
                stop, result, previous = hi, evaluate(hi), evaluate(hi - 1)
                assert result.verdict != 0 and previous.verdict == 0
        rcstop = int(data['rc_stop'][oi])
        assert stop <= rcstop
        rc, _, _ = certificate(*exact_args(data, seen[rcstop]))
        assert int(rc) == int(data['truth'])
        bounds = dict(order=int(oi), stop=stop, verdict=result.verdict)
        if previous is not None:
            bounds['pass_before'] = previous.pass_witness.tolist()
            bounds['fail_before'] = previous.adverse_witness.tolist()
            bounds['adverse_margin'] = previous.adverse_margin
        if result.pass_witness is not None: bounds['pass_at_stop'] = result.pass_witness.tolist()
        boundaries.append(bounds)
        rows.append(dict(**meta, order=int(oi), exact_stop=stop, rc_stop=rcstop,
            exact_checks=int(costs[stop]), rc_checks=int(costs[rcstop]), screen_checks=0 if meta['screen'] else E,
            verdict=result.verdict, truth=int(data['truth']), record_saving=int(costs[rcstop] - costs[stop]),
            exact_fraction=float(costs[stop] / E), rc_fraction=float(costs[rcstop] / E),
            offline_uncached_calls=calls, offline_uncached_seconds=seconds, offline_evaluations=evaluations,
            solver_unknowns=unknowns, earliest_proved=True))
    dest.mkdir(exist_ok=True, parents=True)
    key = meta['table_id']
    pd.DataFrame(rows).to_csv(dest / f'{key}.csv', index=False)
    (dest / f'{key}_traces.json').write_text(json.dumps(traces, separators=(',', ':')))
    (dest / f'{key}_boundary_witnesses.json').write_text(json.dumps(boundaries, separators=(',', ':')))
    elapsed = perf_counter() - start
    print(f'{key}: n={len(rows)}, RC={100*np.mean([r["rc_fraction"] for r in rows]):.4f}%, '
          f'Exact={100*np.mean([r["exact_fraction"] for r in rows]):.4f}%, offline={elapsed:.2f}s', flush=True)
    return rows


def online(meta, order_ids, dest):
    """True prefix-by-prefix replay: no cache, no RC bound, no truth hint."""
    data = load(meta['table_id'])
    N, E, K, ep, bf, fc = map(int, data['public'])
    rows = []
    for oi in order_ids:
        order = data['orders'][oi]
        for method in ['RC-VoA', 'Exact MILP', 'Pre-screen + full']:
            seen = np.zeros((7, 8), dtype=np.int64)
            calls = unknowns = queried = 0
            solve_seconds = 0.0
            tic = perf_counter()
            # Fixed public pre-screen shared by all methods.
            screen = min(K, int((data['cells'] // 4 == 1)[(data['cells'] & 2) == 0].sum())) < bf - fc
            assert screen == meta['screen']
            if screen:
                verdict, prefix = -1, 0
            else:
                for prefix in range(len(order) + 1):
                    args = exact_args(data, seen)
                    ts = perf_counter()
                    if method == 'RC-VoA':
                        verdict, _, _ = certificate(*args)
                        verdict = int(verdict); calls += 1
                    elif method == 'Exact MILP':
                        result = ec.verify(*args)
                        verdict = result.verdict; calls += result.calls; unknowns += result.unknowns
                    elif prefix == len(order):
                        verdict = 1 if np.all(ec.full_margins(seen, ep, bf, fc, K) >= 0) else -1
                        calls += 1
                    else: verdict = 0
                    solve_seconds += perf_counter() - ts
                    if verdict: break
                    if prefix < len(order):
                        addition = data['group_counts'][order[prefix]]
                        seen += addition; queried += int(addition.sum())
                else: raise AssertionError('Full records must resolve')
            elapsed = perf_counter() - tic
            assert verdict == int(data['truth'])
            if method == 'RC-VoA': assert prefix == int(data['rc_stop'][oi])
            rows.append(dict(**meta, order=int(oi), method=method, stop=prefix, checks=queried,
                             verdict=verdict, wall_seconds=elapsed, verifier_seconds=solve_seconds,
                             calls=calls, unknowns=unknowns))
    dest.mkdir(exist_ok=True, parents=True)
    pd.DataFrame(rows).to_csv(dest / f'{meta["table_id"]}.csv', index=False)
    print('ONLINE', meta['table_id'], [(m, round(g.wall_seconds.mean(), 4))
          for m, g in pd.DataFrame(rows).groupby('method')], flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['pilot', 'full', 'online'])
    parser.add_argument('--orders', type=int, default=3, help='Online first N archived orders (default 3)')
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    table = pd.read_csv(HERE / 'inputs/tables.csv')
    if args.mode == 'pilot':
        table = table[~table.screen].sort_values('table_id').groupby(['stage', 'cohort', 'feasible'], as_index=False).head(1)
        ids, dest = [0, 49, 99], HERE / 'pilot'
    elif args.mode == 'full': ids, dest = list(range(100)), HERE / 'offline'
    else:
        table = table[~table.screen]
        ids, dest = list(range(args.orders)), HERE / 'online'
    rows = []
    for meta in table.to_dict('records'):
        saved = dest / f'{meta["table_id"]}.csv'
        if args.resume and saved.exists():
            old = pd.read_csv(saved)
            assert sorted(old.order.unique()) == ids
            rows.extend(old.to_dict('records'))
            continue
        rows.extend((online if args.mode == 'online' else offline)(meta, ids, dest))
    frame = pd.DataFrame(rows)
    frame.to_csv(HERE / f'{args.mode}_results.csv', index=False)
    receipt = dict(mode=args.mode, rows=len(rows), tables=len(table), orders=ids,
                   python=platform.python_version(), scipy=scipy.__version__, numpy=np.__version__,
                   platform=platform.platform(), processor=platform.processor(),
                   complete=True, disagreements=int((frame.verdict != np.where(frame.feasible, 1, -1)).sum()))
    (HERE / f'{args.mode}_receipt.json').write_text(json.dumps(receipt, indent=2))
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == '__main__': main()
