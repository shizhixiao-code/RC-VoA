"""Exact target sensitivity and system/cohort summaries from frozen cell counts.

Run: python analyze_targets.py
Requires NumPy and SciPy (milp). Outputs are deterministic feasibility counts;
solver witnesses may differ when multiple optimal allocations exist.
"""
from pathlib import Path
import csv, hashlib, itertools, json, platform, re, time
from collections import Counter, defaultdict
import numpy as np
import scipy
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.optimize._highspy._core import HIGHS_VERSION_MAJOR, HIGHS_VERSION_MINOR, HIGHS_VERSION_PATCH

HERE=Path(__file__).resolve().parent
OUT=HERE/'results'
RATES={'Melanoma':95,'Invasive SCC':95,'SCC NOS':95,'BCC':90,'SCCIS':90}

def ceil_pct(rate,n): return (int(rate)*int(n)+99)//100
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write_json(path,x): path.write_text(json.dumps(x,indent=2),encoding='utf-8')
def write_csv(path,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def arrays(tab):
    labels=np.array([c[0] for c in tab['cells']])
    b,p,r,n=(np.array([c[j] for c in tab['cells']],dtype=int) for j in range(1,5))
    ca=~np.isin(labels,['Benign','Excluded']);be=labels=='Benign'
    return labels,b,p,r,n,ca,be

def spec(profile,tab):
    s=dict(PROTOCOL['reference']);s.update(profile['change'])
    s['cap']=min(100,int(tab['cap'])+s['absolute_cap_delta_pp'])
    return s

def satisfies(tab,x,s):
    labels,b,p,r,n,ca,be=arrays(tab);N=int(n.sum());final=n*b+(p-b)*x
    if np.any(x<0) or np.any(x>n) or int(x.sum())>s['imaging_pct']*N//100:return False
    if 100*int(final[ca].sum())<s['pooled_pct']*int(n[ca].sum()):return False
    for label,rate in RATES.items():
        use=labels==label
        if 100*int(final[use].sum())<(rate+s['subtype_delta_pp'])*int(n[use].sum()):return False
    for d in [b,r]:
        if int(final[ca&(d==1)].sum())!=int(n[ca&(d==1)].sum()):return False
        if 100*int(final[be].sum())>(100-s['fp_reduction_pct'])*int((n*d)[be].sum()):return False
        if 100*int((n*d).sum()-final.sum())<s['flag_reduction_pp']*N:return False
    return 100*int(final.sum())<=s['cap']*N

def allocation_milp(tab,s):
    labels,b,p,r,n,ca,be=arrays(tab);N=int(n.sum());delta=p-b;f=n*b
    rows=[];lo=[];hi=[]
    def add(coeff,lower=-np.inf,upper=np.inf):rows.append(coeff);lo.append(lower);hi.append(upper)
    add(np.ones(len(n)),upper=s['imaging_pct']*N//100)
    add(delta*ca,lower=ceil_pct(s['pooled_pct'],n[ca].sum())-int(f[ca].sum()))
    for label,rate in RATES.items():
        mask=labels==label
        add(delta*mask,lower=ceil_pct(rate+s['subtype_delta_pp'],n[mask].sum())-int(f[mask].sum()))
    for d in [b,r]:
        mask=ca&(d==1)
        add(delta*mask,lower=int(n[mask].sum()-f[mask].sum()))
        add(delta*be,upper=(100-s['fp_reduction_pct'])*int((n*d)[be].sum())//100-int(f[be].sum()))
        add(delta,upper=int((n*d).sum())-ceil_pct(s['flag_reduction_pp'],N)-int(f.sum()))
    add(delta,upper=s['cap']*N//100-int(f.sum()))
    result=milp(np.ones(len(n)),integrality=np.ones(len(n)),bounds=Bounds(np.zeros(len(n)),n),
                constraints=LinearConstraint(np.array(rows),lo,hi),options={'time_limit':30,'mip_rel_gap':0.0})
    if result.status==2:return -1,None
    if result.status!=0:raise RuntimeError((tab['table_id'],s,result.status,result.message))
    x=np.rint(result.x).astype(int)
    assert np.max(np.abs(result.x-x))<1e-6
    assert satisfies(tab,x,s),(tab['table_id'],s,x)
    return int(x.sum()),x.tolist()

def closed_form(tab,s):
    labels,b,p,r,n,ca,be=arrays(tab);N=int(n.sum());union=b|r
    up=(b==0)&(p==1);down=(b==1)&(p==0)
    w=int(n[ca&(r==1)&(b==0)].sum());v=int(n[ca&(r==1)&(b==0)&(p==0)].sum())
    needs=[];union_needs=[];support=[]
    for label,rate in RATES.items():
        mask=labels==label;target=ceil_pct(rate+s['subtype_delta_pp'],n[mask].sum())
        needs.append(max(0,target-int((n*b)[mask].sum())))
        union_needs.append(max(0,target-int((n*union)[mask].sum())))
        support.append(int(n[mask&up].sum()))
    total=ceil_pct(s['pooled_pct'],n[ca].sum())
    tstar=w+max(sum(union_needs),max(0,total-int((n*union)[ca].sum())))
    F=int((n*b).sum());FR=int((n*r).sum())
    fmax=min(s['cap']*N//100,F-ceil_pct(s['flag_reduction_pp'],N),FR-ceil_pct(s['flag_reduction_pp'],N))
    BB=int((n*b)[be].sum());BR=int((n*r)[be].sum())
    m=BB-min((100-s['fp_reduction_pct'])*BB//100,(100-s['fp_reduction_pct'])*BR//100)
    required=max(m,tstar+F-fmax);cm=int(n[ca&up].sum());bp=int(n[be&down].sum())
    u=int(n[(labels=='Excluded')&down].sum());K=s['imaging_pct']*N//100
    margins=[-v,min(a-b for a,b in zip(support,needs)),cm-tstar,bp-m,bp+u-required,K-tstar-required]
    good=min(margins)>=0
    return {'feasible':good,'minimum_images':tstar+required if good else -1,'t_star':tstar,'C_minus':cm,'B_plus':bp,'m':m,'r':required,'K':K,'margins':margins}

def metadata(tab):
    name=tab['table_id']
    if tab['family']=='source_SP':
        system=name[len(tab['cohort'])+1:];backbone,fusion=system.rsplit('_',1)
    else:
        m=re.fullmatch(r'(?:Repaired_)?(?:MRA|UQ)_(.+)_b[1-5]_(?:frozen|q)(?:30|50)',name)
        assert m,name
        system=m[1];backbone='';fusion=''
    return dict(system=system,backbone=backbone,fusion=fusion)

def exhaustive_validation():
    rng=np.random.default_rng(20261108);checked=feasible=0
    for N in range(2,9):
        for rep in range(6):
            cells=[]
            for _ in range(N):
                b=int(rng.integers(2));r=b if rep%2==0 else int(rng.integers(2))
                cells.append([str(rng.choice(list(RATES)+['Benign','Benign','Excluded'])),b,int(rng.integers(2)),r,1])
            tab={'table_id':f'synthetic_{N}_{rep}','cap':int(rng.choice([30,50,100])),'cells':cells}
            for profile in PROTOCOL['profiles']:
                s=spec(profile,tab)
                literal=min((sum(x) for x in itertools.product([0,1],repeat=N) if satisfies(tab,np.array(x),s)),default=-1)
                derived=closed_form(tab,s);solved,_=allocation_milp(tab,s)
                assert literal==solved==derived['minimum_images'],(tab,s,literal,solved,derived)
                checked+=1;feasible+=literal>=0
    return {'judgments':checked,'feasible':feasible,'seed':20261108}

def main():
    start=time.perf_counter();OUT.mkdir(exist_ok=True)
    tables=json.loads((HERE/'tables.json').read_text());assert len(tables)==236
    rows=[];witnesses=[];reference={}
    for profile in PROTOCOL['profiles']:
        for tab in tables:
            s=spec(profile,tab);cert=closed_form(tab,s);minimum,x=allocation_milp(tab,s)
            assert minimum==cert['minimum_images'],(tab['table_id'],profile,minimum,cert)
            if profile['name']=='reference':
                assert cert['feasible']==tab['expected_feasible'] and minimum==tab['expected_minimum']
                assert cert['margins']==tab['expected_margins']
                reference[tab['table_id']]=cert['feasible']
            else:assert not reference[tab['table_id']] or cert['feasible']
            row=dict(table_id=tab['table_id'],family=tab['family'],cohort=tab['cohort'],batch=tab['batch'],original_cap=tab['cap'],N=tab['N'],**metadata(tab),profile=profile['name'],**s,
                     **{k:v for k,v in cert.items() if k!='margins'},margins=';'.join(map(str,cert['margins'])))
            rows.append(row)
            if x is not None:witnesses.append(dict(table_id=tab['table_id'],profile=profile['name'],selected_cell_counts=x))
        print(profile['name'],{f:sum(r['feasible'] for r in rows if r['profile']==profile['name'] and r['family']==f) for f in ['source_SP','primary','repaired']},flush=True)
    write_csv(OUT/'target_sensitivity_all.csv',rows);write_json(OUT/'allocation_witnesses.json',witnesses)
    summary=[]
    for profile in PROTOCOL['profiles']:
        for family in ['source_SP','primary','repaired']:
            sub=[r for r in rows if r['profile']==profile['name'] and r['family']==family]
            summary.append(dict(profile=profile['name'],family=family,tables=len(sub),feasible=sum(r['feasible'] for r in sub),newly_feasible=sum(r['feasible'] and not reference[r['table_id']] for r in sub)))
    write_csv(OUT/'target_summary.csv',summary)
    refs=[r for r in rows if r['profile']=='reference']
    grouped=defaultdict(list)
    for r in refs:grouped[(r['family'],r['cohort'],r['system'],r['original_cap'])].append(r)
    systems=[]
    for (family,cohort,system,cap),sub in sorted(grouped.items()):
        systems.append(dict(family=family,cohort=cohort,system=system,cap=cap,tables=len(sub),feasible=sum(r['feasible'] for r in sub),feasible_batches=';'.join(str(r['batch']) for r in sub if r['feasible'])))
    write_csv(OUT/'system_cohort_summary.csv',systems)
    strata=[]
    for family in ['primary','repaired']:
        for cohort in ['MRA','UQ']:
            for cap in [30,50]:
                sub=[r for r in refs if r['family']==family and r['cohort']==cohort and r['original_cap']==cap]
                strata.append(dict(family=family,cohort=cohort,cap=cap,tables=len(sub),feasible=sum(r['feasible'] for r in sub)))
    write_csv(OUT/'denominator_reconciliation.csv',strata)
    small=exhaustive_validation()
    audit=dict(status='PASS',tables=len(tables),profiles=len(PROTOCOL['profiles']),independent_milp_judgments=len(rows),
               reference_agreement=len(reference),integer_witnesses_checked=len(witnesses),small_table_enumeration=small,
               relaxations_monotone=True,unknown_or_timeout=0,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
               highs=f'{HIGHS_VERSION_MAJOR}.{HIGHS_VERSION_MINOR}.{HIGHS_VERSION_PATCH}',elapsed_seconds=time.perf_counter()-start,
               hashes={p.name:sha(p) for p in [HERE/'tables.json',HERE/'PROTOCOL.json',Path(__file__)]})
    write_json(OUT/'AUDIT.json',audit)
    print(json.dumps(audit,indent=2))

PROTOCOL=json.loads((HERE/'PROTOCOL.json').read_text())
if __name__=='__main__':main()
