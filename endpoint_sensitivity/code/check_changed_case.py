"""Check the only changed configuration with independent lesion-variable MILPs."""
from pathlib import Path
import csv,json
from collections import Counter
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from run_analysis import altered,SETTINGS,LABELS

ROOT=Path(__file__).resolve().parents[1]
TABLE='Repaired_UQ_localheads_b3_q50'

def main():
    rows=[]
    with np.load(ROOT/'inputs/tables'/f'{TABLE}.npz',allow_pickle=False) as z:
        removed=[{'diagnosis':str(z['diagnosis'][i]),'B':int((z['cells'][i]>>2)&1),'P':int((z['cells'][i]>>1)&1),'R':int(z['cells'][i]&1)} for i in np.flatnonzero(z['labels']==7)]
        removed_counts=Counter((x['diagnosis'],x['B'],x['P'],x['R']) for x in removed)
        with (ROOT/'results/changed_case_removed_cells.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['diagnosis','B','P','R','records'])
            w.writerows([(*key,n) for key,n in sorted(removed_counts.items())])
        for setting in SETTINGS:
            labels,cells,_,_=altered(z,setting)
            n=len(labels);b=((cells&4)!=0).astype(int);p=((cells&2)!=0).astype(int);r=((cells&1)!=0).astype(int)
            delta=p-b;ca=labels<6;nc=labels==6
            for budget in [15,20]:
                coeff=[];lower=[];upper=[]
                def add(a,lo=-np.inf,hi=np.inf):coeff.append(a);lower.append(lo);upper.append(hi)
                add(np.ones(n),hi=budget*n//100)
                add(delta*ca,lo=(90*int(ca.sum())+99)//100-int(b[ca].sum()))
                for j,rate in enumerate([95,95,95,90,90]):
                    sub=labels==j;add(delta*sub,lo=(rate*int(sub.sum())+99)//100-int(b[sub].sum()))
                for baseline in [b,r]:
                    protected=ca&(baseline==1)
                    add(delta*protected,lo=int(protected.sum()-b[protected].sum()))
                    add(delta*nc,hi=90*int(baseline[nc].sum())//100-int(b[nc].sum()))
                    add(delta,hi=int(baseline.sum())-(10*n+99)//100-int(b.sum()))
                add(delta,hi=50*n//100-int(b.sum()))
                model=np.array(coeff)
                result=milp(np.ones(n),integrality=np.ones(n),bounds=Bounds(np.zeros(n),np.ones(n)),constraints=LinearConstraint(model,lower,upper),options={'time_limit':30,'mip_rel_gap':0.0})
                assert result.status in [0,2]
                feasible=result.status==0
                assert feasible==(setting!='endpoint_complete_case')
                if feasible:
                    x=np.rint(result.x).astype(int);assert np.max(np.abs(x-result.x))<1e-6
                    assert np.all(model@x>=np.array(lower)) and np.all(model@x<=np.array(upper))
                    final=b+(p-b)*x
                    assert np.all(final[ca&((b==1)|(r==1))]==1)
                    minimum=int(x.sum());assert minimum==164
                else:minimum=-1
                downgrade=((b==1)&(p==0));t=4
                max_flags=min(50*n//100,int(b.sum())-(10*n+99)//100,int(r.sum())-(10*n+99)//100)
                needed=t+int(b.sum())-max_flags
                rows.append({'table_id':TABLE,'setting':setting,'imaging_pct':budget,'N':n,'K':budget*n//100,'variables':n,'solver_status':result.status,'feasible':feasible,'minimum_images':minimum,'required_rescues':t,'available_total_downgrades':int(downgrade.sum()),'required_total_downgrades':needed,'downgrade_margin':int(downgrade.sum())-needed})
    with (ROOT/'results/changed_case_lesion_milp.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    out={'status':'PASS','independent_lesion_milps':len(rows),'changed_configuration':TABLE,'results':rows}
    (ROOT/'validation/CHANGED_CASE_AUDIT.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
