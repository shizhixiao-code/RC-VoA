"""Conditional paired group bootstrap, with common weights across configurations."""
from pathlib import Path
import csv,hashlib,json,platform,time
from collections import defaultdict
import numpy as np
import scipy
from reference_allocation import closed_form,allocation_milp
HERE=Path(__file__).resolve().parent
LABELS=['Melanoma','Invasive SCC','SCC NOS','BCC','SCCIS','Other cancer','Benign','Excluded']
C=np.arange(8);B=(C&4)!=0;P=(C&2)!=0;R=(C&1)!=0;A=B|R
UP=~B&P;DOWN=B&~P;MAND=~B&R;MISSING=MAND&~P
RATES=np.array([95,95,95,90,90])

def write_csv(name,rows):
    with (HERE/'results'/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def audit(x,cap):
    """x[...,diagnostic category,decision cell], including excluded category 7."""
    N=x.sum(axis=(-2,-1));ca=x[...,:6,:].sum(axis=-2);ben=x[...,6,:]
    nj=x[...,:5,:].sum(axis=-1);req=(RATES*nj+99)//100
    needs=np.maximum(0,req-(x[...,:5,:]*B).sum(axis=-1))
    un=np.maximum(0,req-(x[...,:5,:]*A).sum(axis=-1))
    t=(ca*MAND).sum(axis=-1)+np.maximum(un.sum(axis=-1),np.maximum(0,(90*ca.sum(axis=-1)+99)//100-(ca*A).sum(axis=-1)))
    BB=(ben*B).sum(axis=-1);BR=(ben*R).sum(axis=-1)
    m=np.maximum((BB+9)//10,BB-9*BR//10)
    FB=(x*B).sum(axis=(-2,-1));FR=(x*R).sum(axis=(-2,-1))
    fmax=np.minimum(np.minimum(cap*N//100,FB-(N+9)//10),FR-(N+9)//10)
    rr=np.maximum(m,t+FB-fmax)
    margins=np.stack([-(ca*MISSING).sum(axis=-1),((x[...,:5,:]*UP).sum(axis=-1)-needs).min(axis=-1),(ca*UP).sum(axis=-1)-t,(ben*DOWN).sum(axis=-1)-m,(ben*DOWN).sum(axis=-1)+(x[...,7,:]*DOWN).sum(axis=-1)-rr,N//5-t-rr],axis=-1)
    return (margins>=0).all(axis=-1),margins

def table(x,cap,name):
    return {'table_id':name,'cap':cap,'cells':[[LABELS[j],int(B[c]),int(P[c]),int(R[c]),int(x[j,c])] for j,c in np.argwhere(x>0)]}

def spec(cap):return dict(imaging_pct=20,pooled_pct=90,subtype_delta_pp=0,fp_reduction_pct=10,flag_reduction_pp=10,cap=cap)

def main():
    start=time.perf_counter();protocol=json.loads((HERE/'PROTOCOL.json').read_text());n=protocol['replicates']
    rng=np.random.default_rng(protocol['seed']);meta=json.loads((HERE/'inputs/configurations.json').read_text())
    records=defaultdict(list)
    for row in meta:records[row['stratum']].append(row)
    all_verdicts=[];observed=[];weights={};validation=[];unit_rows=[]
    with np.load(HERE/'inputs/paired_group_counts.npz',allow_pickle=False) as data:
        for key,items in records.items():
            counts=data[key];G=counts.shape[2]
            w=rng.multinomial(G,np.full(G,1/G),size=n);weights[key]=w
            unit_rows.append({'stratum':key,'groups':G,'lesions':int(counts[0,0].sum()),'configurations':10})
            for i,entry in enumerate(items):
                cap=entry['stages'][0]['cap'];original=counts[i].sum(axis=1)
                truth,margins=audit(original,cap)
                assert truth.tolist()==[s['feasible'] for s in entry['stages']]
                sampled=np.einsum('rg,sglc->rslc',w,counts[i],optimize=True)
                decisions,mm=audit(sampled,cap)
                assert np.all(decisions[:,1]>=decisions[:,0])
                for stage in range(2):
                    scalar=closed_form(table(original[stage],cap,entry['stages'][stage]['table_id']),spec(cap))
                    assert scalar['margins']==margins[stage].tolist()
                    observed.append(dict(**entry['stages'][stage]))
                # Two deterministically selected paired replicates per configuration.
                for k in [0,n-1]:
                    for stage in range(2):
                        tab=table(sampled[k,stage],cap,f'{key}_{i}_{k}_{stage}')
                        scalar=closed_form(tab,spec(cap));assert scalar['margins']==mm[k,stage].tolist()
                        minimum,_=allocation_milp(tab,spec(cap))
                        assert minimum==scalar['minimum_images']
                        validation.append({'stratum':key,'config_index':i,'replicate':k,'stage':stage,'minimum_images':minimum})
                all_verdicts.append(decisions)
            print(key,'paired bootstrap completed',flush=True)
    vv=np.stack(all_verdicts,axis=1) # draws, configurations, primary/repaired
    oo=np.array([[x['stages'][j]['feasible'] for j in range(2)] for x in meta],dtype=int)
    assert oo.sum(axis=0).tolist()==[1,14]
    summary=[]
    subsets={'all':np.ones(len(meta),bool)}
    for co in ['MRA','UQ']:
        subsets[co]=np.array([e['stages'][0]['cohort']==co for e in meta])
    for cap in [30,50]:
        subsets[f'cap_{cap}']=np.array([e['stages'][0]['cap']==cap for e in meta])
    for label,mask in subsets.items():
        primary=100*vv[:,mask,0].mean(axis=1);repaired=100*vv[:,mask,1].mean(axis=1);diff=repaired-primary
        lo,hi=np.quantile(diff,[.025,.975],method='linear')
        summary.append({'stratum':label,'configurations':int(mask.sum()),'primary_feasible':int(oo[mask,0].sum()),'repaired_feasible':int(oo[mask,1].sum()),'observed_gain_pp':float(100*(oo[mask,1]-oo[mask,0]).mean()),'bootstrap_mean_gain_pp':float(diff.mean()),'percentile_95_low_pp':float(lo),'percentile_95_high_pp':float(hi),'replicates':n})
    write_csv('paired_gain_summary.csv',summary)
    write_csv('bootstrap_draws.csv',[{'replicate':k,'primary_feasible':int(vv[k,:,0].sum()),'repaired_feasible':int(vv[k,:,1].sum()),'gain_pp':int(vv[k,:,1].sum())-int(vv[k,:,0].sum())} for k in range(n)])
    write_csv('reference_table_verdicts.csv',observed);write_csv('sampling_units.csv',unit_rows);write_csv('allocation_milp_validation.csv',validation)
    strata=defaultdict(list)
    for i,e in enumerate(meta):
        s=e['stages'][0];strata[(s['cohort'],s['mode'],s['cap'])].append((s['batch'],*oo[i]))
    detail=[]
    for (co,mode,cap),vals in sorted(strata.items()):
        vals=sorted(vals)
        detail.append({'cohort':co,'mode':mode,'cap':cap,'batches':5,'primary_passes':sum(int(v[1]) for v in vals),'repaired_passes':sum(int(v[2]) for v in vals),'primary_by_batch':';'.join(str(int(v[1])) for v in vals),'repaired_by_batch':';'.join(str(int(v[2])) for v in vals),'gain_by_batch':';'.join(str(int(v[2]-v[1])) for v in vals)})
    write_csv('configuration_fold_stability.csv',detail)
    np.savez_compressed(HERE/'results/shared_bootstrap_weights.npz',**weights)
    np.savez_compressed(HERE/'results/bootstrap_verdicts.npz',verdicts=vv)
    result={'status':'PASS','python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'replicates':n,'seed':protocol['seed'],'cohort_batch_strata':10,'groups':sum(u['groups'] for u in unit_rows),'original_verdicts_matched':200,'bootstrap_table_verdicts':n*200,'independent_allocation_milps':len(validation),'monotone_paired_gain':True,'overall':summary[0],'elapsed_seconds':time.perf_counter()-start,'hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'PROTOCOL.json',HERE/'inputs/paired_group_counts.npz',HERE/'inputs/configurations.json',Path(__file__)]}}
    (HERE/'results/AUDIT.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

if __name__=='__main__':main()
