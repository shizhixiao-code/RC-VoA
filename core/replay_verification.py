"""Recalculate the reported verification comparisons from archived ledgers."""
import csv,json,statistics
from pathlib import Path
from reproduce import load_tables, full_margins
HERE=Path(__file__).resolve().parent
def replay():
    truth={t['table_id']:full_margins(t)['feasible'] for t in load_tables()}
    with (HERE/'verification_ledgers.csv').open(newline='',encoding='utf-8') as f:rows=list(csv.DictReader(f))
    assert len(rows)==23600
    for r in rows:
        r['feasible']=r['feasible']=='True';r['screen_reject']=r['screen_reject']=='True'
        for key in ['order','checks','eligible','verdict']:r[key]=int(r[key])
        assert r['feasible']==truth[r['table_id']]==(r['verdict']==1)
        assert 0<=r['checks']<=r['eligible']
        r['saved']=100*(1-r['checks']/r['eligible'])
    selectors={
        'All primary':lambda r:r['family']=='primary',
        'Primary queried infeasible':lambda r:r['family']=='primary' and not r['feasible'] and not r['screen_reject'],
        'Primary feasible':lambda r:r['family']=='primary' and r['feasible'],
        'Repaired queried infeasible':lambda r:r['family']=='repaired' and not r['feasible'] and not r['screen_reject'],
        'Repaired feasible':lambda r:r['family']=='repaired' and r['feasible']}
    summary=[]
    for scope,predicate in selectors.items():
        sub=[r for r in rows if predicate(r)]
        counts={k:sum(r['table_id']==k for r in sub) for k in {r['table_id'] for r in sub}}
        assert set(counts.values())=={100}
        saved=statistics.fmean(r['saved'] for r in sub)
        row=dict(scope=scope,tables=len(counts),orders=100,avoided_pct=saved,checked_pct=100-saved)
        if scope=='Repaired feasible':
            row['necessary_query_bound_pct']=statistics.fmean(100*float(r['lower_bound'])/r['eligible'] for r in sub)
        summary.append(row)
    expected=[92.02307135616829,74.15851801094676,0.027149321266967796,67.18604361609476,0.01598810545011112]
    assert all(abs(r['avoided_pct']-v)<1e-9 for r,v in zip(summary,expected))
    report=dict(status='PASS',full_record_agreements=len(rows),summary=summary,
                scope='Archived query-cost ledger replay; exact full-table margins recomputed from independent aggregate inputs')
    (HERE/'results').mkdir(exist_ok=True)
    (HERE/'results/VERIFICATION_REPLAY.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
if __name__=='__main__':replay()
