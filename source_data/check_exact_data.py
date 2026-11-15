"""Recompute Figure 4 / Table SVI summaries from packaged ledgers."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parent/'exact_verification'
d=pd.read_csv(R/'full_results.csv');s=pd.read_csv(R/'summary.csv');t=pd.read_csv(R/'online_results.csv')
assert len(d)==20000 and d.table_id.nunique()==200
assert not d[['table_id','order']].duplicated().any()
assert d.groupby('table_id').size().eq(100).all()
assert d.verdict.eq(d.truth).all() and d.solver_unknowns.eq(0).all()
assert d.exact_checks.le(d.rc_checks).all()
assert d[d.screen].groupby('stage').table_id.nunique().eq(72).all()
assert (d.loc[d.screen,['rc_checks','exact_checks']]==0).all().all()
n=d[~d.screen];assert len(n)==5600
assert len(n[~n.feasible])==4100
assert n[~n.feasible].exact_stop.eq(n[~n.feasible].rc_stop).all()
assert len(t)==1680 and not t.unknowns.any()
for row in s.itertuples():
 p=n[n.stage.eq(row.stage)&n.feasible.eq(row.feasible)]
 order=p.groupby('order')[['rc_fraction','exact_fraction']].mean()*100
 assert p.table_id.nunique()==row.configurations
 for m in ['rc','exact']:
  v=order[m+'_fraction']
  for name,actual in [('mean_pct',v.mean()),('p10',v.quantile(.1)),('p90',v.quantile(.9))]:
   assert np.isclose(actual,getattr(row,m+'_'+name),rtol=0,atol=1e-10)
 for method,m in [('RC-VoA','rc'),('Exact MILP','exact'),('Pre-screen + full','full')]:
  v=t[t.stage.eq(row.stage)&t.feasible.eq(row.feasible)&t.method.eq(method)]
  assert np.isclose(v.wall_seconds.mean(),getattr(row,m+'_online_mean_s'),rtol=0,atol=1e-10)
f=n[n.stage.eq('Repaired')&n.feasible]
assert (f.record_saving>0).sum()==419 and f.record_saving.max()==103
assert np.isclose(f.record_saving.mean(),5.472142857142857)
print('Figure 4 / Table SVI: PASS; 20,000 verdicts, 4,100 identical infeasible stops, 560 online pairs.')
