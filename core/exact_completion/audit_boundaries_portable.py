"""Independent allocation MILPs validate saved opposite-completion witnesses.

This model contains selection counts only, no deficit, ceiling or max auxiliaries.
It applies literal protection, sensitivity, flag, FP and image-cap constraints.
"""
from pathlib import Path
import hashlib
import json
from time import perf_counter
import numpy as np
import pandas as pd
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import csc_array
from exact_completion import B, P, R, RATES, full_margins
from run_experiment import load

HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]


def allocation(counts, ep, bf, fc, K):
    counts = np.asarray(counts, dtype=np.int64)
    ub = np.r_[counts.ravel(), ep]
    rows, lower, upper = [], [], []
    delta = P.astype(np.int64) - B.astype(np.int64)
    def add(row, lo=-np.inf, hi=np.inf):
        rows.append(row); lower.append(lo); upper.append(hi)
    # Identity protection separately within each category / public decision cell.
    for j in range(6):
        for c in range(8):
            if (B[c] or R[c]) and counts[j, c]:
                row = np.zeros(57, dtype=np.int64); row[8*j+c] = delta[c]
                add(row, (1-int(B[c]))*counts[j, c])
    for labs, rate in [(list(range(6)), 90)] + [([j], int(r)) for j, r in enumerate(RATES)]:
        row = np.zeros((7, 8), dtype=np.int64)
        row[labs] = 100 * delta
        current = int((counts[labs] * B).sum())
        total = int(counts[labs].sum())
        add(np.r_[row.ravel(), 0], rate * total - 100 * current)
    row = np.zeros((7, 8), dtype=np.int64); row[6] = 10 * delta
    fpb, fpr = int((counts[6] * B).sum()), int((counts[6] * R).sum())
    add(np.r_[row.ravel(), 0], hi=-fpb)
    add(np.r_[row.ravel(), 0], hi=9*fpr-10*fpb)
    add(np.r_[np.tile(delta, 7), -1], hi=fc-bf)
    add(np.ones(57, dtype=np.int64), hi=K)
    mat = np.array(rows)
    result = milp(np.zeros(57), integrality=np.ones(57), bounds=Bounds(np.zeros(57), ub),
        constraints=LinearConstraint(csc_array(mat.astype(float)), lower, upper),
        options={'time_limit':10, 'mip_rel_gap':0.0})
    if result.status == 2: return False
    assert result.status == 0, result.message
    v = np.rint(result.x).astype(np.int64)
    assert np.max(np.abs(v-result.x)) < 1e-5
    assert np.all(v >= 0) and np.all(v <= ub)
    assert np.all(mat@v >= lower) and np.all(mat@v <= upper)
    return True


def joint_completion_allocation(seen, unresolved, ep, bf, fc, K):
    """Independent existential model: completion x and acquisition z together."""
    seen = np.asarray(seen, dtype=np.int64)
    upper_x = seen + unresolved[None, :]
    lb = np.r_[seen.ravel(), np.zeros(57, dtype=np.int64)]
    ub = np.r_[upper_x.ravel(), upper_x.ravel(), ep]
    rows, lows, highs_ = [], [], []
    delta = P.astype(np.int64) - B.astype(np.int64)
    def add(row, lo=-np.inf, hi=np.inf):
        rows.append(row); lows.append(lo); highs_.append(hi)
    for c in range(8):
        row = np.zeros(113, dtype=np.int64)
        row[np.arange(7)*8+c] = 1
        total = int(seen[:,c].sum()+unresolved[c])
        add(row,total,total)
    for j in range(7):
        for c in range(8):
            i = 8*j+c
            row = np.zeros(113, dtype=np.int64);row[i]=-1;row[56+i]=1
            add(row,hi=0)
            if j < 6 and (B[c] or R[c]):
                row = np.zeros(113,dtype=np.int64)
                row[i]=int(B[c])-1;row[56+i]=delta[c]
                add(row,0)
    for labs,rate in [(list(range(6)),90)]+[([j],int(r)) for j,r in enumerate(RATES)]:
        x = np.zeros((7,8),dtype=np.int64);z=x.copy()
        x[labs]=100*B.astype(np.int64)-rate;z[labs]=100*delta
        add(np.r_[x.ravel(),z.ravel(),0],0)
    for reference in [B,R]:
        x=np.zeros((7,8),dtype=np.int64);z=x.copy()
        x[6]=10*B.astype(np.int64)-9*reference.astype(np.int64);z[6]=10*delta
        add(np.r_[x.ravel(),z.ravel(),0],hi=0)
    add(np.r_[np.zeros(56,dtype=np.int64),np.tile(delta,7),-1],hi=fc-bf)
    add(np.r_[np.zeros(56,dtype=np.int64),np.ones(57,dtype=np.int64)],hi=K)
    mat=np.array(rows)
    result=milp(np.zeros(113),integrality=np.ones(113),bounds=Bounds(lb,ub),
        constraints=LinearConstraint(csc_array(mat.astype(float)),lows,highs_),
        options={'time_limit':10,'mip_rel_gap':0.0})
    if result.status==2:return False
    assert result.status==0,result.message
    v=np.rint(result.x).astype(np.int64)
    assert np.max(np.abs(v-result.x))<1e-5 and np.all(v>=lb) and np.all(v<=ub)
    assert np.all(mat@v>=lows) and np.all(mat@v<=highs_)
    return True


def main():
    table = pd.read_csv(HERE/'inputs/tables.csv')
    total_pairs = unique_witnesses = full_tables = negative_stop_solves = 0
    rows = []
    start = perf_counter()
    for meta in table.to_dict('records'):
        data = load(meta['table_id'])
        N, E, K, ep, bf, fc = map(int, data['public'])
        assert allocation(data['group_counts'].sum(0), ep, bf, fc, K) == meta['feasible']
        full_tables += 1
        path = HERE/'offline'/f'{meta["table_id"]}_boundary_witnesses.json'
        boundaries = json.loads(path.read_text())
        cache = {}
        stop_cache = set()
        count = 0
        for row in boundaries:
            if not meta['feasible']:
                stop_seen = data['group_counts'][data['orders'][row['order'], :row['stop']]].sum(0)
                stop_key = stop_seen.tobytes()
                if stop_key not in stop_cache:
                    assert not joint_completion_allocation(stop_seen,data['cell_totals']-stop_seen.sum(0),ep,bf,fc,K)
                    negative_stop_solves += 1;stop_cache.add(stop_key)
            if row['stop'] == 0: continue
            seen = data['group_counts'][data['orders'][row['order'], :row['stop']-1]].sum(0)
            for field, expected in [('pass_before', True), ('fail_before', False)]:
                x = np.array(row[field], dtype=np.int64)
                assert np.all(x >= seen)
                assert np.array_equal(x.sum(0), data['cell_totals'])
                key = x.tobytes()
                if key not in cache:
                    cache[key] = allocation(x, ep, bf, fc, K)
                    unique_witnesses += 1
                assert cache[key] == expected
                assert bool(np.all(full_margins(x,ep,bf,fc,K)>=0)) == expected
            total_pairs += 1; count += 1
        rows.append(dict(table_id=meta['table_id'], boundary_pairs=count, unique_witnesses=len(cache),
                         independent_negative_stop_models=len(stop_cache)))
    preparation = json.loads((HERE/'PREPARATION.json').read_text())
    for p, h in preparation['input_hashes'].items():
        assert hashlib.sha256((HERE/'inputs'/p).read_bytes()).hexdigest() == h, p
    pd.DataFrame(rows).to_csv(HERE/'portable_boundary_audit_by_table.csv',index=False)
    report = dict(status='PASS', full_tables=full_tables, opposite_completion_pairs=total_pairs,
                  independently_solved_unique_witnesses=unique_witnesses,
                  independently_solved_negative_stop_models=negative_stop_solves,
                  independent_milps=full_tables+unique_witnesses+negative_stop_solves, disagreement=0,
                  frozen_inputs_unchanged=True, wall_seconds=perf_counter()-start)
    (HERE/'PORTABLE_BOUNDARY_AUDIT.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__': main()
