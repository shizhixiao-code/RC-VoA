"""Validate archived target/system summaries and input/code hashes (standard library)."""
from pathlib import Path
from collections import Counter,defaultdict
import csv,hashlib,json

def check(root):
    root=Path(root);results=root/'results'
    def read(name):
        with (results/name).open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))
    audit=json.loads((results/'AUDIT.json').read_text())
    for name,digest in audit['hashes'].items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,name
    assert audit['status']=='PASS' and audit['reference_agreement']==236
    assert audit['independent_milp_judgments']==1652 and audit['small_table_enumeration']['judgments']==294
    tables={t['table_id']:t for t in json.loads((root/'tables.json').read_text())}
    rows=read('target_sensitivity_all.csv')
    assert len(rows)==1652 and len({(r['table_id'],r['profile']) for r in rows})==1652
    grouped=defaultdict(list);by_system=defaultdict(list);by_cap=defaultdict(list)
    for r in rows:
        assert r['table_id'] in tables
        grouped[(r['profile'],r['family'])].append(r)
        if r['profile']=='reference':
            t=tables[r['table_id']]
            assert (r['feasible']=='True')==t['expected_feasible']
            assert int(r['minimum_images'])==t['expected_minimum']
            by_system[(r['family'],r['cohort'],r['system'],r['original_cap'])].append(r)
            by_cap[(r['family'],r['cohort'],r['original_cap'])].append(r)
    for r in read('target_summary.csv'):
        sub=grouped[r['profile'],r['family']]
        assert len(sub)==int(r['tables']) and sum(s['feasible']=='True' for s in sub)==int(r['feasible'])
        assert sum(s['feasible']=='True' and not tables[s['table_id']]['expected_feasible'] for s in sub)==int(r['newly_feasible'])
    for r in read('system_cohort_summary.csv'):
        sub=by_system[r['family'],r['cohort'],r['system'],r['cap']]
        assert len(sub)==int(r['tables']) and sum(s['feasible']=='True' for s in sub)==int(r['feasible'])
    for r in read('denominator_reconciliation.csv'):
        sub=by_cap[r['family'],r['cohort'],r['cap']]
        assert len(sub)==int(r['tables']) and sum(s['feasible']=='True' for s in sub)==int(r['feasible'])
    assert sum(r['feasible']=='True' for r in rows if r['profile']=='reference' and r['family']=='repaired')==14
    assert sum(r['feasible']=='True' for r in rows if r['profile']=='flag_reduction_5' and r['family']=='repaired')==26
    batches=read('reference_batch_summary.csv')
    assert len(batches)==5
    for b in batches:
        sub=[r for r in rows if r['profile']=='reference' and r['family']=='repaired' and r['batch']==b['batch']]
        assert len(sub)==int(b['configurations'])==20
        assert sum(r['feasible']=='True' for r in sub)==int(b['feasible'])
    assert [int(b['feasible']) for b in batches]==[7,0,4,3,0]
    return {'status':'PASS','table_profile_rows':len(rows),'system_cohort_rows':len(by_system),'hashes_checked':len(audit['hashes']),'reference_batch_counts':[int(b['feasible']) for b in batches]}

if __name__=='__main__':
    root=Path(__file__).resolve().parent/'target_sensitivity'
    if not root.exists():root=Path(__file__).resolve().parents[1]/'analysis'
    print(json.dumps(check(root),indent=2))
