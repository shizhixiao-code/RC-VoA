"""Validate the restored S5 source data and regenerate the current Table SV rows."""
from pathlib import Path
import csv,json,math
HERE=Path(__file__).resolve().parent

def read(name):
    with (HERE/name).open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))

def number(value):return float(value) if value.strip() else None

def main():
    full=read('Table_SIV_full.csv');op=read('operating_points.csv');cal=read('calibration_thresholds.csv')
    assert len(full)==36 and len(cal)==60
    index={(x['backbone'],x['cohort'],x['setting']):x for x in full};assert len(index)==36
    metrics=['N','excluded','benign_N','flags','FP']+[k for k in full[0] if k.endswith(('_N','_TP'))]
    for row in full:
        if row['setting']=='oracle_cohort':continue
        sub=[x for x in op if x['backbone']==row['backbone'] and x['cohort']==row['cohort'] and x['setting']==row['setting'] and int(x['evaluation_fold'])>=0]
        assert len(sub)==5
        for k in metrics:assert sum(int(x[k]) for x in sub)==int(row[k]),(row['backbone'],row['cohort'],k)
        if row['setting']=='S_local':
            cut=[float(x['threshold']) for x in cal if x['backbone']==row['backbone'] and x['cohort']==row['cohort']]
            assert len(cut)==5
            assert abs(min(cut)-float(row['threshold_min']))<1e-12 and abs(max(cut)-float(row['threshold_max']))<1e-12
    curve=read('frontier.csv');assert len(curve)==27516
    uncertainty=read('uncertainty.csv');assert len(uncertainty)==1728
    for row in uncertainty:
        assert int(row['resamples'])==1000
        assert int(row['valid_resamples'])+int(row['undefined_resamples'])==1000
        lo,hi=number(row['lower']),number(row['upper'])
        if lo is not None and hi is not None:assert lo<=hi
    tables=read('complete_criterion_tables.csv');assert len(tables)==1296
    for row in read('audit_summary.csv'):
        sub=[x for x in tables if all(x[k]==row[k] for k in ['scope','setting','budget'])]
        assert len(sub)==int(row['tables']) and sum(x['feasible']=='True' for x in sub)==int(row['feasible'])
    for budget,want in [('15',96),('20',99)]:
        sub=[x for x in tables if x['scope']=='heldout_fold' and x['setting']=='S_local' and x['budget']==budget]
        assert len(sub)==180 and sum(x['feasible']=='True' for x in sub)==want
    rows=[]
    for co in ['MRA','HIBA','UQ','MCR-SL']:
        for model in ['efficientnet_b0','convnext_tiny','deit3_base_patch16_224']:
            b=index[model,co,'source_S'];f=index[model,co,'S_local'];o=index[model,co,'oracle_cohort']
            rows.append({'cohort':co,'backbone':model,'baseline':'source_S','baseline_FP':int(b['FP']),'final_FP':int(f['FP']),'Delta_FP':int(b['FP'])-int(f['FP']),'source_flags_pct':b['flag_rate'],'oracle_flags_pct':o['flag_rate'],'local_flags_pct':f['flag_rate'],'pooled_TP':f['All cancer_TP'],'pooled_N':f['All cancer_N'],'melanoma_TP':f['Melanoma_TP'],'melanoma_N':f['Melanoma_N'],'all_floors':f['all_observed_floors_met']})
    assert [x['Delta_FP'] for x in rows]==[16,49,45,23,30,34,1576,2216,1599,42,44,44]
    with (HERE/'Table_SV_current.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    report={'status':'PASS','full_table_rows':len(full),'fold_and_pooled_operating_rows':len(op),'calibration_threshold_rows':len(cal),'frontier_rows':len(curve),'bootstrap_interval_rows':len(uncertainty),'complete_audit_rows':len(tables),'Table_SV_baseline_minus_final_rows':len(rows),'legacy_filename_note':'Table_SIV_full.csv is the original filename; its current manuscript display is Table SV.'}
    (HERE/'PORTABLE_VALIDATION.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))

if __name__=='__main__':main()
