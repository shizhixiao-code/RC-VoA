"""Unmodified scalar allocation functions captured from target-sensitivity analysis."""
import numpy as np
from scipy.optimize import Bounds,LinearConstraint,milp
RATES={'Melanoma': 95, 'Invasive SCC': 95, 'SCC NOS': 95, 'BCC': 90, 'SCCIS': 90}

def ceil_pct(rate,n): return (int(rate)*int(n)+99)//100

def arrays(tab):
    labels=np.array([c[0] for c in tab['cells']])
    b,p,r,n=(np.array([c[j] for c in tab['cells']],dtype=int) for j in range(1,5))
    ca=~np.isin(labels,['Benign','Excluded']);be=labels=='Benign'
    return labels,b,p,r,n,ca,be

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