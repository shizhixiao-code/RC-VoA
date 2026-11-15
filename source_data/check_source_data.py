"""Rebuild the manuscript's selector results from frozen aggregate counts."""
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from frozen_audit import ORDER, arrays, full_margins, passes

ROOT = Path(__file__).resolve().parent
DATA = ROOT/'source_data'

def read(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))

def cell_key(row):
    return (row['diagnosis'], int(row['B']), int(row['P']), int(row['R']))

def equal(value, expected, field):
    if value != expected:
        raise ValueError(f'{field}: calculated {value!r}, expected {expected!r}')

def evaluate(table, x):
    _, b, p, r, n, cancer, benign = arrays(table)
    final = n*b + (p-b)*x
    protected = cancer & ((b|r)==1)
    return dict(images=int(x.sum()), quota=int(x[b==0].sum()),
                flags=int(final.sum()), TP=int(final[cancer].sum()), FP=int(final[benign].sum()),
                rescued=int(final[cancer & (b==0)].sum()), lost=int((n-final)[protected].sum()),
                current_lost=int((n-final)[cancer & (b==1)].sum()),
                FP_reduction=int((n*b-final)[benign].sum()), passes=bool(passes(table,x,ORDER)))

def quota_range(table, certificate, volume):
    if not certificate['feasible']:
        return []
    labels,b,p,r,n,cancer,benign = arrays(table)
    neutral=int(n[(b==0)&(p==0)].sum())
    excluded_up=int(n[(b==0)&(p==1)&(labels=='Excluded')].sum())
    excluded_down=int(n[(b==1)&(p==0)&(labels=='Excluded')].sum())
    safe_pos=int(n[(b==1)&((p==1)|~cancer)].sum())
    base_flags=int((n*b).sum())
    max_flags=min(int(n.sum())//2,min(base_flags,int((n*r).sum()))-(int(n.sum())+9)//10)
    valid=[]
    for k in range(certificate['t_star'], min(volume,int(n[b==0].sum()))+1):
        pos=volume-k
        upgrades=max(certificate['t_star'],k-neutral)
        new_fp=max(0,upgrades-certificate['C_minus']-excluded_up)
        ruleouts=min(pos,certificate['B_plus'])
        downgrades=min(pos,certificate['B_plus']+excluded_down)
        if pos<=safe_pos and ruleouts-new_fp>=certificate['m'] and downgrades-upgrades>=base_flags-max_flags:
            valid.append(k)
    return valid

def main():
    outcomes=read(DATA/'all50_configurations_three_policies.csv')
    expected={(r['table_id'],r['policy']):r for r in outcomes}
    equal(len(expected),150,'distinct configuration and policy pairs')
    equal(len(outcomes),len(expected),'unique rows')
    totals=defaultdict(Counter)
    chosen=defaultdict(Counter)
    for row in read(ROOT/'policy_cells.csv'):
        key=(row['table_id'],row['policy'])
        cell=cell_key(row)
        count=int(row['count'])
        totals[key][cell]+=count
        chosen[key][cell]+=count*int(row['selected'])
    equal(set(totals),set(expected),'policy cell coverage')
    tables={}
    observed={}
    stats=defaultdict(Counter)
    for (tid,policy),counts in totals.items():
        row=expected[tid,policy]
        cells=sorted(counts)
        table=dict(table_id=tid, cap=50, cells=[[*c,counts[c]] for c in cells])
        if tid in tables:
            equal(table,tables[tid],'identical frozen inputs across policies')
        tables[tid]=table
        equal(sum(counts.values()),int(row['N']),'N')
        volume=sum(counts.values())//5
        x=np.array([chosen[tid,policy][c] for c in cells])
        result=evaluate(table,x)
        for field,value in result.items():
            equal(value, row[field]=='True' if field=='passes' else int(row[field]), field)
        equal(result['images'],volume,'20% volume')
        cert=full_margins(table)
        valid=quota_range(table,cert,volume)
        equal(bool(valid),row['oracle_exact_feasible']=='True','exact-volume feasibility')
        if valid:
            equal(valid,[int(k) for k in row['quota_values'].split(';')],'complete feasible quotas')
            equal((min(valid),max(valid)),(int(float(row['quota_min'])),int(float(row['quota_max']))),'quota endpoints')
        category=('passed' if result['passes'] else 'support' if min(cert['margins'][:-1])<0
                  else 'capacity' if not cert['feasible'] else 'exact_volume' if not valid
                  else 'quota' if result['quota'] not in valid else 'selection')
        equal(category,row['failure'],'first failed audit level')
        observed[tid,policy]=(result,cert,valid,category)
        stats[policy].update(configurations=1,oracle_feasible=int(bool(valid)),passes=int(result['passes']),
                             conditional_passes=int(bool(valid) and result['passes']))
    for row in read(DATA/'policy_summary.csv'):
        equal(dict(stats[row['policy']]),{k:int(v) for k,v in row.items() if k!='policy'},'policy summary')

    diagnostics=read(DATA/'diagnostic_all10.csv')
    wanted={tid for (tid,pol),(_,_,valid,_) in observed.items() if pol=='Transition-aware' and valid}
    equal({r['table_id'] for r in diagnostics},wanted,'Fig. 3 configuration IDs')
    equal(len(diagnostics),len(wanted),'Fig. 3 unique rows')
    failures=Counter()
    for row in diagnostics:
        result,cert,valid,category=observed[row['table_id'],'Transition-aware']
        for field in ['quota','rescued']:
            equal(result[field],int(row[field]),field)
        equal((min(valid),max(valid)),(int(float(row['quota_min'])),int(float(row['quota_max']))),'Fig. 3 interval')
        equal(cert['t_star'],int(row['t_star']),'Fig. 3 rescue target')
        equal(cert['C_minus']-cert['t_star'],int(row['rescue_slack']),'Fig. 3 slack')
        equal(category,row['failure'],'Fig. 3 failure')
        failures[category]+=1

    witness=defaultdict(Counter)
    for row in read(ROOT/'witness_cells.csv'):
        witness[row['table_id']][cell_key(row)]+=int(row['count'])
    witnesses=read(DATA/'tableIII_three_witnesses.csv')
    equal(set(witness),{row['table_id'] for row in witnesses},'Table III witnesses')
    for row in witnesses:
        tid=row['table_id']
        table=tables[tid]
        result,cert,valid,_=observed[tid,'Transition-aware']
        cells=[tuple(c[:4]) for c in table['cells']]
        equal(set(witness[tid])-set(cells),set(),'witness cells in source table')
        oracle=evaluate(table,np.array([witness[tid][c] for c in cells]))
        equal(oracle['passes'],True,'minimum witness passes')
        equal(oracle['images'],cert['minimum_images'],'minimum acquisitions')
        for field,column in [('images','minimum_K'),('rescued','oracle_rescued'),('lost','oracle_protected_misses'),('FP_reduction','oracle_FP_reduction')]:
            equal(oracle[field],int(row[column]),column)
        for field,column in [('images','K'),('rescued','policy_rescued'),('lost','policy_protected_misses'),('current_lost','policy_B_cancellations'),('FP_reduction','policy_FP_reduction'),('flags','policy_flags'),('quota','quota')]:
            equal(result[field],int(row[column]),column)
    return dict(status='PASS',configurations=len(tables),matched_policy_rows=len(observed),
                policies={p:dict(s) for p,s in stats.items()},figure3_rows=len(diagnostics),
                first_failures=dict(failures),tableIII_minimum_witnesses=len(witnesses))

if __name__=='__main__':
    print(json.dumps(main(),indent=2))


if __name__ == "__main__":
    import runpy
    from pathlib import Path
    runpy.run_path(str(Path(__file__).with_name("check_exact_data.py")), run_name="__main__")


if __name__ == "__main__":
    from check_target_data import check
    print(json.dumps(check(Path(__file__).resolve().parent / "target_sensitivity"), indent=2))
