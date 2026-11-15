"""Functions from the accompanying reproduce.py."""
import numpy as np

RATES = {'Melanoma':95, 'Invasive SCC':95, 'SCC NOS':95, 'BCC':90, 'SCCIS':90}

ORDER = ['imaging','pooled','subtypes','protect_B','protect_R','FP','relative_flags','absolute_flags']

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
