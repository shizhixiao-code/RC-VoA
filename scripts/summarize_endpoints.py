"""Reconstruct current SI VI-B and VIII from original endpoint records/results."""
from pathlib import Path
from collections import Counter
import argparse,csv,json

def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(p,rows):
    with p.open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data-dir',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True)
    a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
    inv=read(a.data_dir/'results/raw_diagnosis_inventory.csv')
    mapping={(r['cohort'],r['provider_diagnosis']):int(r['reference_endpoint']) for r in inv}
    counts=[];bases=[]
    for co in ['MRA','HIBA','UQ','MCR-SL']:
        rows=read(a.data_dir/f'evidence/raw_field_slices/{co}__source_fields.csv')
        field={'MRA':'label_status','HIBA':'diagnosis_confirm_type','UQ':'diagnosis_basis','MCR-SL':'histopathology'}[co]
        n=Counter()
        for r in rows:
            diagnosis=r.get('raw_pathology',r.get('diagnosis'));y=mapping[co,diagnosis]
            n[('cancer' if y==1 else 'noncancer' if y==0 else 'excluded',r[field])]+=1
        ak=sum(int(r['records']) for r in inv if r['cohort']==co and r['AK_alternative_affected']=='True')
        excluded=sum(int(r['records']) for r in inv if r['cohort']==co and int(r['reference_endpoint'])<0)
        counts.append(dict(cohort=co,N=len(rows),explicit_AK=ak,other_excluded=excluded-ak,complete_case_N=len(rows)-excluded))
        bases.extend(dict(cohort=co,endpoint=k[0],field=field,recorded_basis=k[1],n=v) for k,v in sorted(n.items()))
    assert sum(r['N'] for r in counts)==9168
    assert sum(r['explicit_AK'] for r in counts)==92
    assert sum(r['other_excluded'] for r in counts)==107
    assert sum(r['complete_case_N'] for r in counts)==8969
    rows=read(a.data_dir/'results/feasibility_summary.csv');table=[]
    for setting,expected in [('reference',14),('AK_noncancer',14),('endpoint_complete_case',13)]:
        for budget in [15,20]:
            found={r['family']:int(r['feasible']) for r in rows if r['setting']==setting and int(r['imaging_pct'])==budget}
            assert found=={'source_SP':27,'primary':1,'repaired':expected}
            table.append(dict(setting=setting,imaging_pct=budget,**found))
    write(a.output_dir/'Table_SVI_B.csv',table)
    write(a.output_dir/'Table_SVIII_counts.csv',counts)
    write(a.output_dir/'Table_SVIII_recorded_bases.csv',bases)
    receipt={'status':'PASS','scope':'Archived-result and raw-field summarization; use endpoints suite for allocation/MILP replay','table_SVI_B':table,'external_counts':counts}
    (a.output_dir/'CHECK.json').write_text(json.dumps(receipt,indent=2),encoding='utf8')
    print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
