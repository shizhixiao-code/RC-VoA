"""Replay two predeclared endpoint/exclusion settings on 236 frozen tables."""
from pathlib import Path
from collections import Counter, defaultdict
import csv,hashlib,json,platform,time
import numpy as np
import scipy
from reference_allocation import closed_form,allocation_milp,satisfies

ROOT=Path(__file__).resolve().parents[1]
LABELS=['Melanoma','Invasive SCC','SCC NOS','BCC','SCCIS','Other cancer','Benign','Excluded']
AK={'malignant- ak','actinic keratosis','ak'}
SETTINGS=['reference','AK_noncancer','endpoint_complete_case']
MARGINS=['mandatory_rescue','subtype_rescue','total_rescue','FP_downgrade','flag_downgrade','imaging_capacity']

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(name,rows):
    with (ROOT/'results'/name).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def altered(z,setting):
    labels=z['labels'].copy();keep=np.ones(len(labels),bool)
    ak=np.array([str(d).lower() in AK for d in z['diagnosis']])
    assert np.all(labels[ak]==7)
    if setting=='AK_noncancer':labels[ak]=6
    if setting=='endpoint_complete_case':keep=labels!=7
    if setting=='AK_noncancer':assert np.array_equal(labels[~ak],z['labels'][~ak])
    assert np.array_equal(labels<6,z['labels']<6)
    return labels[keep],z['cells'][keep],z['group'][keep],keep
def counts(labels,cells):
    x=np.zeros((8,8),dtype=np.int64);np.add.at(x,(labels,cells),1)
    return x
def tabular(x,entry):
    return {'table_id':entry['table_id'],'cap':entry['cap'],'cells':[[LABELS[j],int((c>>2)&1),int((c>>1)&1),int(c&1),int(x[j,c])] for j,c in np.argwhere(x>0)]}
def stats(labels,cells,groups):
    b=(cells&4)>0;p=(cells&2)>0;r=(cells&1)>0
    return dict(N=len(labels),groups=len(set(groups)),known=int((labels<7).sum()),cancer=int((labels<6).sum()),noncancer=int((labels==6).sum()),excluded=int((labels==7).sum()),
                melanoma=int((labels==0).sum()),invasive_SCC=int((labels==1).sum()),SCC_NOS=int((labels==2).sum()),BCC=int((labels==3).sum()),SCCIS=int((labels==4).sum()),other_cancer=int((labels==5).sum()),
                C_minus=int(((labels<6)&~b&p).sum()),B_plus=int(((labels==6)&b&~p).sum()),excluded_downgrades=int(((labels==7)&b&~p).sum()),
                baseline_flags=int(b.sum()),reference_flags=int(r.sum()),baseline_FP=int(((labels==6)&b).sum()),reference_FP=int(((labels==6)&r).sum()))

def main():
    start=time.perf_counter()
    plan=json.loads((ROOT/'protocol/PLAN.json').read_text())
    assert [s['id'] for s in plan['settings']]==SETTINGS
    input_audit=json.loads((ROOT/'validation/INPUT_AUDIT.json').read_text())
    for rel,h in input_audit['portable_input_hashes'].items():assert sha(ROOT/rel)==h,rel
    entries=json.loads((ROOT/'inputs/tables.json').read_text())
    detail=[];denominators={};witnesses=[];milp_count=0
    for entry in entries:
        with np.load(ROOT/'inputs/tables'/entry['file'],allow_pickle=False) as z:
            assert hashlib.sha256(z['cells'].tobytes()).hexdigest()==entry['decision_sha256']
            for setting in SETTINGS:
                labels,cells,groups,keep=altered(z,setting)
                assert np.array_equal(cells,z['cells'][keep])
                x=counts(labels,cells);table=tabular(x,entry);sums=stats(labels,cells,groups)
                if entry['family']=='source_SP':
                    dkey=(setting,entry['cohort'])
                    fields={k:sums[k] for k in ['N','groups','known','cancer','noncancer','excluded','melanoma','invasive_SCC','SCC_NOS','BCC','SCCIS','other_cancer']}
                    if dkey in denominators:assert denominators[dkey]==fields
                    denominators[dkey]=fields
                budgets=[10,15,20] if entry['family']=='source_SP' else [15,20]
                for budget in budgets:
                    spec=dict(imaging_pct=budget,pooled_pct=90,subtype_delta_pp=0,fp_reduction_pct=10,flag_reduction_pp=10,cap=entry['cap'])
                    cf=closed_form(table,spec);minimum,witness=allocation_milp(table,spec)
                    assert minimum==cf['minimum_images'],(entry['table_id'],setting,budget)
                    assert (minimum>=0)==cf['feasible'];milp_count+=1
                    if witness is not None:
                        assert satisfies(table,np.array(witness),spec)
                        witnesses.append({'table_id':entry['table_id'],'setting':setting,'imaging_pct':budget,'cells':table['cells'],'selected_by_cell':witness,'minimum_images':minimum})
                    margin=cf.pop('margins');cf.pop('C_minus');cf.pop('B_plus')
                    failure='feasible' if cf['feasible'] else 'support_failure' if min(margin[:5])<0 else 'capacity_only'
                    row={k:entry[k] for k in ['table_id','family','cohort','mode','batch','cap']}
                    row.update(setting=setting,imaging_pct=budget,**sums,**cf,required_images=cf['t_star']+cf['r'],failure=failure)
                    row.update({name:int(v) for name,v in zip(MARGINS,margin)})
                    detail.append(row)
    write('table_results.csv',detail)
    write('cohort_denominators.csv',[dict(setting=s,cohort=c,**v) for (s,c),v in denominators.items()])
    (ROOT/'results/allocation_witnesses.json').write_text(json.dumps(witnesses,indent=2),encoding='utf-8')
    ref={(x['table_id'],x['imaging_pct']):x for x in detail if x['setting']=='reference'}
    changes=[]
    for x in detail:
        if x['setting']=='reference':continue
        base=ref[x['table_id'],x['imaging_pct']]
        assert x['C_minus']==base['C_minus'] and x['t_star']==base['t_star']
        changes.append({**{k:x[k] for k in ['table_id','family','cohort','mode','batch','cap','setting','imaging_pct']},
                        'reference_feasible':base['feasible'],'alternative_feasible':x['feasible'],
                        'changed_verdict':base['feasible']!=x['feasible'],'reference_failure':base['failure'],'alternative_failure':x['failure'],
                        **{f'reference_{k}':base[k] for k in ['N','K','required_images','B_plus','excluded_downgrades','m','r']},
                        **{f'alternative_{k}':x[k] for k in ['N','K','required_images','B_plus','excluded_downgrades','m','r']}})
    write('paired_table_changes.csv',changes)
    summary=[];failure_summary=[];stratified=[]
    for setting in SETTINGS:
        for family in ['source_SP','primary','repaired']:
            for budget in ([10,15,20] if family=='source_SP' else [15,20]):
                part=[x for x in detail if (x['setting'],x['family'],x['imaging_pct'])==(setting,family,budget)]
                expected_n=36 if family=='source_SP' else 100;assert len(part)==expected_n
                feasible=sum(x['feasible'] for x in part)
                if setting=='reference':assert feasible==(0 if budget==10 else {'source_SP':27,'primary':1,'repaired':14}[family])
                changeset=[x for x in changes if (x['setting'],x['family'],x['imaging_pct'])==(setting,family,budget)]
                summary.append({'setting':setting,'family':family,'imaging_pct':budget,'tables':len(part),'feasible':feasible,'feasible_pct':100*feasible/len(part),
                                'newly_feasible':sum(x['alternative_feasible'] and not x['reference_feasible'] for x in changeset),'newly_infeasible':sum(x['reference_feasible'] and not x['alternative_feasible'] for x in changeset),
                                'support_failure':sum(x['failure']=='support_failure' for x in part),'capacity_only':sum(x['failure']=='capacity_only' for x in part)})
                failure_summary.append({'setting':setting,'family':family,'imaging_pct':budget,'tables':len(part),**{k:sum(x[k]<0 for x in part) for k in MARGINS}})
                for co in sorted({x['cohort'] for x in part}):
                    for cap in sorted({x['cap'] for x in part if x['cohort']==co}):
                        sub=[x for x in part if x['cohort']==co and x['cap']==cap]
                        stratified.append({'setting':setting,'family':family,'imaging_pct':budget,'cohort':co,'cap':cap,'tables':len(sub),'feasible':sum(x['feasible'] for x in sub)})
    write('feasibility_summary.csv',summary);write('failure_margins_nonexclusive.csv',failure_summary);write('cohort_cap_summary.csv',stratified)
    paired=[]
    for setting in SETTINGS:
        for budget in [15,20]:
            p=next(x for x in summary if (x['setting'],x['family'],x['imaging_pct'])==(setting,'primary',budget))
            q=next(x for x in summary if (x['setting'],x['family'],x['imaging_pct'])==(setting,'repaired',budget))
            paired.append({'setting':setting,'imaging_pct':budget,'configurations':100,'primary_feasible':p['feasible'],'repaired_feasible':q['feasible'],'observed_gain_pp':q['feasible']-p['feasible']})
    write('repair_gain.csv',paired)
    preservation=json.loads((ROOT/'validation/PRESERVED_CURRENT_VERSION.json').read_text())
    old=Path(preservation['original_root'])
    if old.exists():
        for name,h in preservation['files'].items():assert sha(old/name)==h,name
    audit={'status':'PASS','settings':3,'alternative_settings':2,'reference_tables':236,'table_setting_budget_checks':len(detail),'independent_allocation_milps':milp_count,'integer_witnesses_checked':len(witnesses),
           'all_reference_counts_reproduced':True,'all_original_cancer_labels_and_rescue_support_unchanged':True,'model_or_threshold_refits':0,'manuscript_modified':False,'protected_files_unchanged':len(preservation['files']) if old.exists() else 'Original project unavailable; included frozen inputs hash-checked',
           'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'elapsed_seconds':time.perf_counter()-start,'plan_sha256':sha(ROOT/'protocol/PLAN.json'),'script_sha256':sha(Path(__file__)),
           'inputs_sha256':{p.relative_to(ROOT).as_posix():sha(p) for p in (ROOT/'inputs').rglob('*') if p.is_file()}}
    assert len(detail)==milp_count==1524
    (ROOT/'validation/ANALYSIS_AUDIT.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps({'status':'PASS','milps':milp_count,'main_20pct_results':[x for x in summary if x['imaging_pct']==20],'repair_gain':paired},indent=2),flush=True)

if __name__=='__main__':main()
