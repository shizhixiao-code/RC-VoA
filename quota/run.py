"""Exploratory fixed-ranking quota intervention, with archived provenance checks."""
from pathlib import Path
import json,hashlib,shutil,sys
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
WORK=ROOT.parent
# Packaged replay only; no original workstation or recovery paths are accessed.
EVAL=ROOT/'inputs'
PREP=ROOT/'recovered/prepared'
SEALFILE=ROOT/'inputs/SELECTION_SEAL.json'

def copy(src,dst):
    if src.resolve()!=dst.resolve():shutil.copy2(src,dst)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for d in ['inputs','results','recovered/prepared']: (ROOT/d).mkdir(parents=True,exist_ok=True)
plan=dict(status='POST_HOC_EXPLORATORY',scope='All 50 repaired configurations at 20% imaging / 50% flag cap; all three frozen rankings',
    intervention='Keep within-pool ordering, fixed K, decisions and endpoints. Sweep every integer quota allowed by pool sizes. Also project original quota to the archived oracle feasible interval.',
    interpretation='Outcome-informed diagnostic upper bound within frozen rankings, not a deployable policy or new independent validation.',
    no_retraining=True,no_new_selector=True)
planfile=ROOT/'PROTOCOL.json'
if not planfile.exists():planfile.write_text(json.dumps(plan,indent=2),encoding='utf-8')
seal=json.loads(SEALFILE.read_text())
done=json.loads((EVAL/'EVALUATION_COMPLETE.json').read_text())
provenance=[]
for f in ['case_results.csv.gz','configuration_results.csv','EVALUATION_COMPLETE.json']:
    p=EVAL/f
    if f in done['files']:assert sha(p)==done['files'][f]
    copy(p,ROOT/'inputs'/f);provenance.append(dict(source=str(p),copy=str(ROOT/'inputs'/f),sha256=sha(p)))
copy(SEALFILE,ROOT/'inputs/SELECTION_SEAL.json')
cfg=pd.read_csv(ROOT/'inputs/configuration_results.csv')
cfg=cfg[(cfg.family=='repaired')&(cfg.audit_flag_cap==50)]
assert len(cfg)==150 and cfg.table_id.nunique()==50
cases=pd.read_csv(ROOT/'inputs/case_results.csv.gz',keep_default_na=False,low_memory=False)
cases=cases[cases.table_id.isin(cfg.table_id)&(cases.cap==50)]
for name in cfg.table_id.unique():
    src=PREP/name;dr=json.loads((src/'DONE.json').read_text());assert sha(src/'DONE.json')==seal['receipts'][name]
    dest=ROOT/'recovered/prepared'/name;dest.mkdir(exist_ok=True)
    for f in ['DONE.json','probabilities.npz','selection.npz','policy.json']:
        p=src/f
        if f in dr['files']:assert sha(p)==dr['files'][f]
        copy(p,dest/f)
        provenance.append(dict(source=str(p),copy=str(dest/f),sha256=sha(p)))
    with np.load(src/'probabilities.npz',allow_pickle=False) as pr,np.load(src/'selection.npz',allow_pickle=False) as ss:
        for pol in ['Transition-aware','Proximity','Clinical-risk']:
            cr=cases[(cases.table_id==name)&(cases.policy==pol)].set_index('id').loc[pr['id']]
            assert np.array_equal(cr['rank'],pr[pol]) and np.array_equal(cr.selected,ss[pol])
rates={'All cancer':90,'Melanoma':95,'Invasive SCC':95,'SCC NOS':95,'BCC':90,'SCCIS':90}
sweep=[];summaries=[]
for (name,pol),frame in cases.groupby(['table_id','policy'],sort=True):
    row=cfg[(cfg.table_id==name)&(cfg.policy==pol)].iloc[0]
    frame=frame.reset_index(drop=True);n=len(frame);K=int(row.K)
    b=frame.B.to_numpy(bool);p=frame.P.to_numpy(bool);r=frame.R.to_numpy(bool);y=frame.y.to_numpy();c=y==1;be=y==0
    masks={t:(c if t=='All cancer' else c&(frame.subtype.to_numpy()==t)) for t in rates}
    required={t:(rate*int(masks[t].sum())+99)//100 for t,rate in rates.items()}
    flagcap=min(n//2,int(b.sum())-(n+9)//10,int(r.sum())-(n+9)//10)
    fpcap=9*min(int((b&be).sum()),int((r&be).sum()))//10
    def evaluate(s):
        final=np.where(s,p,b)
        deficits={t:max(0,required[t]-int((final&m).sum())) for t,m in masks.items()}
        lost=int((c&(b|r)&~final).sum());fp=int((be&final).sum());flags=int(final.sum())
        clinical=lost+max(deficits['All cancer'],sum(v for t,v in deficits.items() if t!='All cancer'))
        return dict(TP=int((c&final).sum()),FP=fp,flags=flags,rescued=int((c&~b&final).sum()),lost=lost,
            current_lost=int((c&b&~final).sum()),clinical_gap=clinical,flag_gap=max(0,flags-flagcap),FP_gap=max(0,fp-fpcap),
            passes=clinical==0 and flags<=flagcap and fp<=fpcap)
    original=evaluate(frame.selected.to_numpy(bool))
    for k,v in original.items():assert v==row[k],(name,pol,k,v,row[k])
    neg=frame[frame.B==0].sort_values(['rank','id']).index.to_numpy();pos=frame[frame.B==1].sort_values(['rank','id']).index.to_numpy()
    actual=int(row.quota);s=np.zeros(n,bool);s[neg[:actual]]=True;s[pos[:K-actual]]=True
    assert np.array_equal(s,frame.selected.to_numpy(bool))
    local=[]
    for k in range(max(0,K-len(pos)),min(K,len(neg))+1):
        s=np.zeros(n,bool);s[neg[:k]]=True;s[pos[:K-k]]=True
        assert s.sum()==K
        m=evaluate(s);rec=dict(table_id=name,policy=pol,quota=k,**m);local.append(rec);sweep.append(rec)
    good=[v['quota'] for v in local if v['passes']]
    best=min(local,key=lambda v:(not v['passes'],v['clinical_gap'],v['flag_gap']+v['FP_gap'],v['FP'],v['flags'],abs(v['quota']-actual),v['quota']))
    feasible=bool(row.oracle_exact_feasible)
    projected=None
    if feasible:
        lo,hi=int(row.quota_min),int(row.quota_max);projected=min(max(actual,lo),hi)
        assert all(lo<=q<=hi for q in good)
    proj=next((v for v in local if v['quota']==projected),None)
    summaries.append(dict(table_id=name,policy=pol,oracle_feasible=feasible,original_quota=actual,original_pass=original['passes'],
        projected_quota=projected,projected_pass=proj['passes'] if proj else None,
        any_quota_pass=bool(good),passing_quotas=';'.join(map(str,good)),swept_quotas=len(local),
        original_lost=original['lost'],best_quota=best['quota'],best_lost=best['lost'],best_rescued=best['rescued'],
        best_clinical_gap=best['clinical_gap'],best_flag_gap=best['flag_gap'],best_FP_gap=best['FP_gap']))
    if len(summaries)%30==0:print('completed',len(summaries),'/150',flush=True)
pd.DataFrame(sweep).to_csv(ROOT/'results/all_quota_outcomes.csv',index=False)
df=pd.DataFrame(summaries);df.to_csv(ROOT/'results/configuration_summary.csv',index=False)
summary=[]
for pol,g in df.groupby('policy',sort=False):
    f=g[g.oracle_feasible]
    summary.append(dict(policy=pol,configurations=len(g),oracle_feasible=len(f),original_pass=int(g.original_pass.sum()),
        projected_pass=int(f.projected_pass.sum()),any_quota_pass=int(g.any_quota_pass.sum()),feasible_any_quota_pass=int(f.any_quota_pass.sum())))
pd.DataFrame(summary).to_csv(ROOT/'results/summary.csv',index=False)
report=dict(status='PASS',original_metrics_reproduced=150,sealed_ranking_selection_contexts=50,quota_evaluations=len(sweep),summary=summary,provenance=provenance)
(ROOT/'VALIDATION.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='provenance'},indent=2))
