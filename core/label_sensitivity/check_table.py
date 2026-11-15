"""Reconstruct Supplementary Table SVII from the completed P1 summaries (stdlib)."""
from pathlib import Path
import csv,json,hashlib

ROOT=Path(__file__).resolve().parent
def read(name):
    with (ROOT/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def get(rows,family,budget,column):
    found=[r for r in rows if r['setting']=='reference' and r['family']==family and int(r['imaging_pct'])==budget]
    assert len(found)==1
    return int(found[0][column])
def main():
    result=[]
    for name,file,col in [('Original','original_feasibility_summary.csv','feasible'),('AK/KA, updated R','pipeline_feasibility_summary.csv','expanded_source_feasible'),('AK/KA, fixed R','frozen_feasibility_summary.csv','expanded_source_feasible')]:
        rr=read(file)
        primary=[get(rr,'primary',b,col) for b in [15,20]]
        repaired=[get(rr,'repaired',b,col) for b in [15,20]]
        assert primary==[1,1] and repaired[0]==repaired[1]
        result.append(dict(setting=name,source_SP_15=get(rr,'source_SP',15,col),source_SP_20=get(rr,'source_SP',20,col),source_SP_denominator=36,primary=1,repaired=repaired[0],lowflag_denominator=100))
    assert [(r['source_SP_15'],r['source_SP_20'],r['repaired']) for r in result]==[(27,27,14),(24,24,7),(21,24,7)]
    plan=json.loads((ROOT/'protocol/PLAN.json').read_text())['source_alternative']
    assert (plan['reference_N'],plan['N'],plan['added_AK_negative'],plan['added_KA_positive'])==(4717,4978,99,162)
    fit=json.loads((ROOT/'validation/FITTING_AUDIT.json').read_text())
    assert fit['image_network_training_jobs']==118 and fit['fusion_outer_pipelines']==15
    tex=ROOT.parent/'manuscript/supp_table_inclusion_sensitivity.tex'
    if tex.exists():
        body=tex.read_text()
        for r in result:
            label=r['setting'].replace(' R',' $R$')
            assert f"{label} & {r['source_SP_15']}/36 & {r['source_SP_20']}/36 & {r['repaired']}/100" in body
    with (ROOT/'Table_SVII.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(result[0]));w.writeheader();w.writerows(result)
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.csv')}
    (ROOT/'TABLE_AUDIT.json').write_text(json.dumps(dict(status='PASS',setting='reference',budgets=[15,20],clinical_targets='Original reference targets',rows=result,input_and_output_csv_sha256=hashes),indent=2)+'\n')
    print(json.dumps(dict(status='TABLE_SVII_PASS',rows=result)))
if __name__=='__main__':main()
