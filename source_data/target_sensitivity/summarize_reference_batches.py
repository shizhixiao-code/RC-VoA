"""Reproduce the supplementary S4 batch counts from reference-target rows."""
from pathlib import Path
import csv,json
HERE=Path(__file__).resolve().parent
with (HERE/'results/target_sensitivity_all.csv').open(newline='',encoding='utf-8') as f:
    data=list(csv.DictReader(f))
rows=[]
for batch in range(1,6):
    sub=[r for r in data if r['profile']=='reference' and r['family']=='repaired' and int(r['batch'])==batch]
    assert len(sub)==20 and {r['cohort'] for r in sub}=={'MRA','UQ'}
    assert {int(r['original_cap']) for r in sub}=={30,50}
    rows.append(dict(batch=batch,feasible=sum(r['feasible']=='True' for r in sub),configurations=len(sub)))
assert [r['feasible'] for r in rows]==[7,0,4,3,0]
assert sum(r['feasible'] for r in rows)==14
with (HERE/'results/reference_batch_summary.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps({'status':'PASS','reference_repaired_by_batch':rows},indent=2))
