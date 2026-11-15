"""Independent cumulative-count verification of every directly enumerated allocation."""
from pathlib import Path
import pandas as pd
import numpy as np
import json
R=Path(__file__).resolve().parent
c=pd.read_csv(R/'inputs/case_results.csv.gz',keep_default_na=False,low_memory=False)
q=pd.read_csv(R/'results/all_quota_outcomes.csv'); checks=0
rates={'Melanoma':95,'Invasive SCC':95,'SCC NOS':95,'BCC':90,'SCCIS':90}
for (name,policy),alloc in q.groupby(['table_id','policy']):
    x=c[(c.table_id==name)&(c.policy==policy)&(c.cap==50)].copy()
    n=len(x);K=n//5;b=x.B.to_numpy();p=x.P.to_numpy();r=x.R.to_numpy();y=x.y.to_numpy();ca=(y==1);be=(y==0)
    ni=np.flatnonzero(b==0);pi=np.flatnonzero(b==1)
    rank=x['rank'].to_numpy();ids=x.id.to_numpy();ni=ni[np.lexsort((ids[ni],rank[ni]))];pi=pi[np.lexsort((ids[pi],rank[pi]))]
    def count(mask):
        change=(p-b)*mask
        return int((b*mask).sum())+np.r_[0,np.cumsum(change[ni])][alloc.quota.to_numpy()]+np.r_[0,np.cumsum(change[pi])][K-alloc.quota.to_numpy()]
    tp=count(ca);fp=count(be);flags=count(np.ones(n,int));protected=ca&((b+r)>0);lost=int(protected.sum())-count(protected)
    floors=(100*tp>=90*int(ca.sum()))
    for t,rate in rates.items():
        mask=ca&(x.subtype.to_numpy()==t);floors&=100*count(mask)>=rate*int(mask.sum())
    passes=floors&(lost==0)&(10*fp<=9*min(int((b*be).sum()),int((r*be).sum())))&(flags<=min(n//2,int(b.sum())-(n+9)//10,int(r.sum())-(n+9)//10))
    for field,v in [('TP',tp),('FP',fp),('flags',flags),('lost',lost),('passes',passes)]:assert np.array_equal(v,alloc[field].to_numpy()),(name,policy,field)
    checks+=len(alloc)
out={'status':'PASS','independent_prefix_checks':checks,'fields':['TP','FP','flags','lost','passes'],'method':'Cumulative transition changes and rational sensitivity comparisons; independent of direct selection-loop evaluator'}
(R/'results/INDEPENDENT_CHECK.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(out)
