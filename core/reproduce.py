"""Replay exact frozen-table feasibility and post hoc constraint ablations.

Only NumPy and SciPy are needed. No images, model weights, patient IDs or
workspace-specific modules are required. Run --example or --all --check.
"""
from pathlib import Path
import argparse, csv, hashlib, itertools, json, platform
import numpy as np
import scipy
from scipy.optimize import Bounds, LinearConstraint, milp

HERE = Path(__file__).resolve().parent
RATES = {'Melanoma':95, 'Invasive SCC':95, 'SCC NOS':95, 'BCC':90, 'SCCIS':90}
ORDER = ['imaging','pooled','subtypes','protect_B','protect_R','FP','relative_flags','absolute_flags']
ALTERNATIVE = ['imaging','protect_B','protect_R','pooled','subtypes','FP','relative_flags','absolute_flags']

def write_csv(path, rows):
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def load_tables():
    with (HERE/'tables.json').open(encoding='utf-8') as f: tables=json.load(f)
    assert len(tables)==236
    return tables

def arrays(tab):
    cells=tab['cells'];labels=np.array([c[0] for c in cells])
    b,p,r,n=(np.array([c[j] for c in cells],dtype=int) for j in range(1,5))
    ca=~np.isin(labels,['Benign','Excluded']);be=labels=='Benign'
    return labels,b,p,r,n,ca,be

def passes(tab, x, active):
    labels,b,p,r,n,ca,be=arrays(tab);final=n*b+(p-b)*x
    N=int(n.sum());K=N*20//100
    if np.any(x<0) or np.any(x>n):return False
    if 'imaging' in active and x.sum()>K:return False
    if 'pooled' in active and 100*final[ca].sum()<90*n[ca].sum():return False
    if 'subtypes' in active:
        for label,rate in RATES.items():
            use=labels==label
            if 100*final[use].sum()<rate*n[use].sum():return False
    for key,decision in [('protect_B',b),('protect_R',r)]:
        if key in active and final[ca & (decision==1)].sum()!=n[ca & (decision==1)].sum():return False
    if 'FP' in active and any(10*final[be].sum()>9*(n*decision)[be].sum() for decision in [b,r]):return False
    if 'relative_flags' in active and any(10*((n*decision).sum()-final.sum())<N for decision in [b,r]):return False
    if 'absolute_flags' in active and 100*final.sum()>tab['cap']*N:return False
    return True

def solve(tab,active):
    labels,b,p,r,n,ca,be=arrays(tab);change=p-b;N=int(n.sum())
    rows=[];lower=[];upper=[]
    def add(coeff,lo=-np.inf,hi=np.inf):
        rows.append(coeff);lower.append(lo);upper.append(hi)
    def floor(mask,rate):
        add(change*mask,lo=(rate*int(n[mask].sum())+99)//100-int((n*b)[mask].sum()))
    if 'imaging' in active:add(np.ones(len(n)),hi=N*20//100)
    if 'pooled' in active:floor(ca,90)
    if 'subtypes' in active:
        for label,rate in RATES.items():floor(labels==label,rate)
    for key,decision in [('protect_B',b),('protect_R',r)]:
        if key in active:
            mask=ca & (decision==1)
            add(change*mask,lo=int((n*(1-b))[mask].sum()))
    if 'FP' in active:
        cap=min(9*int((n*b)[be].sum())//10,9*int((n*r)[be].sum())//10)
        add(change*be,hi=cap-int((n*b)[be].sum()))
    if 'relative_flags' in active:
        cap=min(int((n*b).sum()),int((n*r).sum()))-(N+9)//10
        add(change,hi=cap-int((n*b).sum()))
    if 'absolute_flags' in active:add(change,hi=tab['cap']*N//100-int((n*b).sum()))
    constraints=LinearConstraint(np.array(rows),lower,upper) if rows else None
    result=milp(np.ones(len(n)),integrality=np.ones(len(n)),bounds=Bounds(np.zeros(len(n)),n),
                constraints=constraints,options={'time_limit':30,'mip_rel_gap':0.0})
    assert result.status in [0,2],(tab['table_id'],result.message)
    if result.status==2:return -1,None
    x=np.rint(result.x).astype(int)
    assert np.max(np.abs(result.x-x))<1e-6 and passes(tab,x,active)
    return int(x.sum()),x

def full_margins(tab):
    labels,b,p,r,n,ca,be=arrays(tab);up=(b==0)&(p==1);down=(b==1)&(p==0)
    union=(b|r);w=int(n[ca&(r==1)&(b==0)].sum())
    v=int(n[ca&(r==1)&(b==0)&(p==0)].sum())
    needs=[];needs_union=[];support=[]
    for label,rate in RATES.items():
        mask=labels==label;target=(rate*int(n[mask].sum())+99)//100
        needs.append(max(0,target-int((n*b)[mask].sum())))
        needs_union.append(max(0,target-int((n*union)[mask].sum())))
        support.append(int(n[mask&up].sum()))
    total=(90*int(n[ca].sum())+99)//100
    t=w+max(sum(needs_union),max(0,total-int((n*union)[ca].sum())))
    N=int(n.sum());F=int((n*b).sum())
    fcap=min(tab['cap']*N//100,F-(N+9)//10,int((n*r).sum())-(N+9)//10)
    m=int((n*b)[be].sum())-min(9*int((n*b)[be].sum())//10,9*int((n*r)[be].sum())//10)
    req=max(m,t+F-fcap);cm=int(n[ca&up].sum());bp=int(n[be&down].sum())
    u=int(n[(labels=='Excluded')&down].sum())
    margins=[-v,min(np.array(support)-needs),cm-t,bp-m,bp+u-req,N//5-t-req]
    return dict(t_star=t,m=m,r=req,C_minus=cm,B_plus=bp,minimum_images=t+req if min(margins[:-1])>=0 else -1,
                feasible=bool(min(margins)>=0),margins=[int(x) for x in margins])

def exhaustive_check():
    rng=np.random.default_rng(20261104);checks=0
    patterns=[set(ORDER[:i]) for i in range(1,9)]+[set(ORDER)-{s} for s in ORDER]
    for N in range(1,8):
        for rep in range(4):
            tab={'table_id':f'synthetic_{N}_{rep}','cap':int(rng.choice([30,50,100])),
                 'cells':[[str(rng.choice(list(RATES)+['Other cancer','Benign','Excluded'])),
                           int(rng.integers(2)),int(rng.integers(2)),int(rng.integers(2)),1] for _ in range(N)]}
            choices=[np.array(x) for x in itertools.product([0,1],repeat=N)]
            for active in patterns:
                exact=min((int(x.sum()) for x in choices if passes(tab,x,active)),default=-1)
                computed,_=solve(tab,active);assert exact==computed,(tab,active,exact,computed)
                checks+=1
    return checks

def run_all(check):
    out=HERE/'results';out.mkdir(exist_ok=True)
    exhaustive=exhaustive_check() if check else None
    tables=load_tables();rows=[];witnesses=[];margins=[]
    for tab in tables:
        cache={}
        def evaluate(active):
            key=tuple(sorted(active))
            if key not in cache:cache[key]=solve(tab,set(active))
            return cache[key]
        for path,order in [('cumulative',ORDER),('protection_first',ALTERNATIVE)]:
            previous=True
            for i in range(1,len(order)+1):
                minimum,x=evaluate(order[:i]);feasible=minimum>=0
                assert previous or not feasible;previous=feasible
                rows.append(dict(table_id=tab['table_id'],family=tab['family'],cohort=tab['cohort'],batch=tab['batch'],cap=tab['cap'],
                    analysis=path,step=i,constraint=order[i-1],active=';'.join(order[:i]),feasible=feasible,minimum_images=minimum))
                if x is not None and path=='cumulative':
                    witnesses.append(dict(table_id=tab['table_id'],step=i,selected_cell_counts=x.tolist()))
        for removed in ORDER:
            active=set(ORDER)-{removed};minimum,_=evaluate(active)
            rows.append(dict(table_id=tab['table_id'],family=tab['family'],cohort=tab['cohort'],batch=tab['batch'],cap=tab['cap'],
                analysis='leave_one_out',step=0,constraint=removed,active=';'.join(s for s in ORDER if s in active),feasible=minimum>=0,minimum_images=minimum))
        full=full_margins(tab);minimum,_=evaluate(ORDER)
        assert full['feasible']==(minimum>=0)==tab['expected_feasible']
        if minimum>=0:assert minimum==full['minimum_images']==tab['expected_minimum']
        assert full['margins']==tab['expected_margins']
        margins.append(dict(table_id=tab['table_id'],family=tab['family'],**full))
    write_csv(out/'all_configurations.csv',rows)
    summary=[]
    keys=sorted(set((r['analysis'],r['step'],r['constraint'],r['family']) for r in rows))
    for analysis,step,constraint,family in keys:
        sub=[r for r in rows if (r['analysis'],r['step'],r['constraint'],r['family'])==(analysis,step,constraint,family)]
        summary.append(dict(analysis=analysis,step=step,constraint=constraint,family=family,tables=len(sub),feasible=sum(r['feasible'] for r in sub)))
    write_csv(out/'summary.csv',summary)
    (out/'cell_witnesses.json').write_text(json.dumps(witnesses,separators=(',',':')))
    (out/'complete_margins.json').write_text(json.dumps(margins,indent=2))
    primary=[m for m in margins if m['family']=='primary']
    failures=[sum(m['margins'][i]<0 for m in primary) for i in range(6)]
    assert failures==[99,95,99,3,79,50],failures
    report=dict(status='PASS',tables=len(tables),configuration_stage_judgments=len(rows),
        exhaustive_small_table_checks=exhaustive,full_criterion_agreement=len(margins),
        primary_nonexclusive_failures=failures,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        data_sha256=hashlib.sha256((HERE/'tables.json').read_bytes()).hexdigest(),
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'AUDIT.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
    for row in summary:
        if row['analysis']=='cumulative':print(row)

def example():
    for tab in load_tables():
        if tab['table_id'] not in ['UQ_multiscale_b3_frozen50','Repaired_UQ_headonly_b1_q50']:continue
        minimum,x=solve(tab,set(ORDER));cert=full_margins(tab)
        assert minimum>=0 and cert['feasible'] and minimum==cert['minimum_images']
        labels,b,p,r,n,ca,be=arrays(tab);final=n*b+(p-b)*x
        print(json.dumps(dict(table_id=tab['table_id'],minimum_images=minimum,
            retained=int(final[ca&((b|r)==1)].sum()),rescued=int(final[ca&(b==0)].sum()),
            false_positive_reduction=int((n*b)[be].sum()-final[be].sum()),
            flags=int(final.sum()),**{k:cert[k] for k in ['t_star','m','r','C_minus','B_plus','margins']}),indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--all',action='store_true');parser.add_argument('--check',action='store_true');parser.add_argument('--example',action='store_true')
    args=parser.parse_args()
    if args.all:run_all(args.check)
    else:example()
