"""Capture exactly Figure 4's frozen decisions, grouped orders and RC costs."""
from pathlib import Path
import ast
import hashlib
import json
import numpy as np
import pandas as pd
from exact_completion import full_margins, RATES, B, P, R, DOWN

HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]
SOURCE = WORK / 'rcvoa_threshold_transport_20261102'
EXT = SOURCE / 'analysis/rescue_verification_extension_20261026'
LABELS = ['Melanoma', 'Invasive SCC', 'SCC NOS', 'BCC', 'SCCIS']


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    dest = HERE / 'inputs'
    dest.mkdir(exist_ok=True)
    hashes = {}
    def record(p):
        hashes[str(p.relative_to(WORK))] = sha(p)
        return p
    # Extract only the original interval functions, without legacy path mutation.
    src = record(SOURCE / 'code/dual_verification.py')
    parsed = ast.parse(src.read_text(encoding='utf-8'))
    functions = [ast.get_source_segment(src.read_text(encoding='utf-8'), n)
                 for n in parsed.body if isinstance(n, ast.FunctionDef) and n.name in ('floor_bounds', 'certificate')]
    header = ('# Unmodified functions captured from archived dual_verification.py.\n'
              'import numpy as np\nfrom exact_completion import RATES as RATE, B as BASE, '
              'R as REF, A as UNION, UP, DOWN, MAND as MANDATORY, MISSING\n\n')
    (HERE / 'rc_interval.py').write_text(header + '\n\n'.join(functions) + '\n', encoding='utf-8')
    from rc_interval import certificate
    primary = pd.read_csv(record(SOURCE / 'results/verification_increment/paired_orders.csv'))
    repaired = pd.read_csv(record(EXT / 'verification/all_10000_paired_orders.csv'))
    figure_ledger = pd.read_csv(record(WORK / 'rcvoa_figures_integrated_20261103/figure_sources/rcvoa_fig4_verification_20261103/source_data/paired_order_records.csv'))
    assert np.array_equal(figure_ledger.checks, repaired.checks)
    assert np.array_equal(figure_ledger.old_checks, repaired.old_checks)
    fits = json.loads(record(SOURCE / 'inputs/convnext_tiny.json').read_text())
    records = []
    for co in ['MRA', 'UQ']:
        with np.load(record(SOURCE / f'inputs/convnext_tiny__{co}.npz'), allow_pickle=False) as z:
            raw = {k: z[k] for k in z.files}
        all_labels = np.where(raw['y'] < 0, -1, np.where(raw['y'] == 0, 6, 5))
        for j, name in enumerate(LABELS):
            all_labels[(raw['y'] == 1) & (raw['fine_subtype'] == name)] = j
        for mode in ['headonly', 'laststage', 'multiscale', 'localheads', 'localfeatures']:
            with np.load(record(SOURCE / f'inputs/{co}__{mode}__planning_decisions.npz'), allow_pickle=False) as z:
                plan = {k: z[k] for k in z.files}
            assert np.array_equal(plan['id'], raw['id'])
            for fold in range(5):
                use = plan['fold'] == fold
                groups, inv = np.unique(raw['group'][use].astype(str), return_inverse=True)
                assert not (set(groups) & set(raw['group'][~use].astype(str)))
                labels = all_labels[use]
                eligible = labels >= 0
                N, E = len(labels), int(eligible.sum())
                K = N // 5
                for cap in [30, 50]:
                    q = cap / 100
                    with np.load(record(SOURCE / f'results/dual_verification/prefixes/{co}__{mode}__fold{fold}__q{q}.npz'), allow_pickle=False) as z:
                        orders = z['orders'].copy()
                        assert np.array_equal(groups, z['groups'].astype(str))
                    assert orders.shape == (100, len(groups))
                    assert np.all(np.sort(orders, axis=1) == np.arange(len(groups)))
                    for stage in ['Primary', 'Repaired']:
                        key = f'{stage}_{co}_{mode}_b{fold+1}_q{cap}'
                        if stage == 'Primary':
                            b, p = plan[f'{q}_base'][use], plan[f'{q}_pair'][use]
                            r = raw['p1'][use] >= fits['S']['p1']
                            expected = primary[(primary.cohort == co) & (primary['mode'] == mode) &
                                               (primary.fold == fold) & (primary.q == q)].sort_values('order')
                        else:
                            with np.load(record(EXT / 'prepared' / f'{key}.npz'), allow_pickle=False) as z:
                                assert np.array_equal(z['indices'], np.flatnonzero(use))
                                assert np.array_equal(z['id'], raw['id'][use])
                                b, p, r = z['base'].copy(), z['pair'].copy(), z['ref'].copy()
                            expected = repaired[repaired.table_id == key].sort_values('order')
                        cells = 4 * b.astype(int) + 2 * p.astype(int) + r.astype(int)
                        counts = np.zeros((len(groups), 7, 8), dtype=np.int64)
                        np.add.at(counts, (inv[eligible], labels[eligible], cells[eligible]), 1)
                        cell_totals = np.bincount(cells[eligible], minlength=8)
                        ep = int(((~eligible) & DOWN[cells]).sum())
                        bf = int(b.sum())
                        fc = min(cap * N // 100, bf - (N + 9) // 10, int(r.sum()) - (N + 9) // 10)
                        screen = min(K, int((b & ~p).sum())) < bf - fc
                        actual = full_margins(counts.sum(0), ep, bf, fc, K)
                        truth = 1 if np.all(actual >= 0) else -1
                        seen = np.concatenate([np.zeros((100, 1, 7, 8), dtype=np.int64), counts[orders].cumsum(1)], axis=1)
                        unresolved = cell_totals - seen.sum(2)
                        verdict, _, _ = certificate(seen, unresolved, ep, bf, fc, K)
                        stop = (verdict != 0).argmax(1)
                        cost = seen.sum((2, 3))
                        assert len(expected) == 100
                        assert np.array_equal(expected.stop_prefix, stop)
                        assert np.array_equal(expected.checks, cost[np.arange(100), stop])
                        assert np.all(expected.verdict == truth)
                        assert np.all(verdict[verdict != 0] == truth)
                        assert np.all((stop == 0) == screen)
                        np.savez_compressed(dest / f'{key}.npz', group_counts=counts, orders=orders,
                            cell_totals=cell_totals, rc_stop=stop, rc_checks=cost[np.arange(100), stop],
                            public=np.array([N, E, K, ep, bf, fc]), truth=np.array(truth),
                            labels=labels, cells=cells, inv=inv)
                        records.append(dict(table_id=key, stage=stage, cohort=co, mode=mode, batch=fold+1,
                            cap=cap, N=N, eligible=E, groups=len(groups), feasible=truth == 1, screen=bool(screen)))
    df = pd.DataFrame(records)
    df.to_csv(dest / 'tables.csv', index=False)
    assert len(df) == 200
    assert df.groupby('stage').screen.sum().to_dict() == {'Primary': 72, 'Repaired': 72}
    assert df.groupby('stage').feasible.sum().to_dict() == {'Primary': 1, 'Repaired': 14}
    manuscript_hashes = {str(p.relative_to(WORK)): sha(p) for root in [WORK/'rcvoa_integrated_copy_20261105/manuscript',
                          WORK/'rcvoa_submission_final_20261106'] for p in root.glob('*') if p.is_file()}
    receipt = dict(source='Exact Figure 4 inputs; no refitting or new outcome information',
                   source_hashes=hashes, input_hashes={p.name: sha(p) for p in dest.glob('*')},
                   manuscript_hashes=manuscript_hashes, tables=200, matched_archived_costs=20000,
                   interval_source_sha256=sha(src), solver_inputs='Only observed counts, unresolved cell sizes, public constants')
    (HERE / 'PREPARATION.json').write_text(json.dumps(receipt, indent=2), encoding='utf-8')
    print(df.groupby(['stage', 'screen', 'feasible']).size().to_string(), flush=True)


if __name__ == '__main__': main()
