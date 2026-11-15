"""Aggregate archived lesion records into paired group/cell sufficient counts."""
from pathlib import Path
import ast,hashlib,io,json,zipfile
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
LABELS=['Melanoma','Invasive SCC','SCC NOS','BCC','SCCIS','Other cancer','Benign','Excluded']

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

with zipfile.ZipFile(ROOT/'original_version/RCVoA_Reproducibility.zip') as z:
    groups={};meta=[];source_hashes={}
    for cohort in ['MRA','UQ']:
        for batch in range(1,6):
            key=f'{cohort}_b{batch}';tensor=[];first=None
            for mode in ['headonly','laststage','multiscale','localheads','localfeatures']:
                for cap in [30,50]:
                    pair=[];info=[]
                    for stage in ['Primary','Repaired']:
                        name=f'{stage}_{cohort}_{mode}_b{batch}_q{cap}'
                        content=z.read(f'RCVoA_Reproducibility/exact_completion/inputs/{name}.npz')
                        source_hashes[name]=hashlib.sha256(content).hexdigest()
                        with np.load(io.BytesIO(content),allow_pickle=False) as a:
                            labels=np.where(a['labels']<0,7,a['labels']);inv=a['inv'];cells=a['cells']
                            if first is None:first=(labels.copy(),inv.copy())
                            else:assert np.array_equal(first[0],labels) and np.array_equal(first[1],inv),name
                            G=int(inv.max())+1
                            counts=np.zeros((G,8,8),dtype=np.int64)
                            np.add.at(counts,(inv,labels,cells),1)
                            assert np.array_equal(counts[:,:7],a['group_counts'])
                            assert counts.sum()==int(a['public'][0])
                            pair.append(counts);info.append({'table_id':name,'stage':stage,'cohort':cohort,'batch':batch,'mode':mode,'cap':cap,'N':int(a['public'][0]),'groups':G,'feasible':bool(int(a['truth'])==1)})
                    # Fixed labels/group identity and baseline/reference decisions are preserved by repair.
                    assert np.array_equal(pair[0].sum(axis=2),pair[1].sum(axis=2))
                    for mask in [(np.arange(8)&4)!=0,(np.arange(8)&1)!=0]:
                        assert np.array_equal(pair[0][:,:,mask].sum(2),pair[1][:,:,mask].sum(2))
                    tensor.append(np.stack(pair));meta.append({'stratum':key,'configuration_index':len(tensor)-1,'stages':info})
            groups[key]=np.stack(tensor) # configuration, stage, group, label, cell
    np.savez_compressed(HERE/'inputs/paired_group_counts.npz',**groups)
(HERE/'inputs/configurations.json').write_text(json.dumps(meta,indent=2))
references=json.loads((ROOT/'analysis/tables.json').read_text())
for record in meta:
    x=groups[record['stratum']][record['configuration_index']].sum(axis=1)
    for j,stage in enumerate(record['stages']):
        co,mode,b,cap=stage['cohort'],stage['mode'],stage['batch'],stage['cap']
        old_id=f'{co}_{mode}_b{b}_frozen{cap}' if j==0 else f'Repaired_{co}_{mode}_b{b}_q{cap}'
        tab=next(t for t in references if t['table_id']==old_id)
        mapped=np.zeros((8,8),dtype=np.int64)
        for label,B,P,R,n in tab['cells']:mapped[LABELS.index(label),4*B+2*P+R]+=n
        assert np.array_equal(mapped,x[j]),stage['table_id']
source=ROOT/'analysis/analyze_targets.py'
tree=ast.parse(source.read_text(encoding='utf-8'))
keep={'ceil_pct','arrays','satisfies','allocation_milp','closed_form'}
code='"""Unmodified scalar allocation functions captured from target-sensitivity analysis."""\nimport numpy as np\nfrom scipy.optimize import Bounds,LinearConstraint,milp\nRATES='+repr({'Melanoma':95,'Invasive SCC':95,'SCC NOS':95,'BCC':90,'SCCIS':90})+'\n\n'
code+='\n\n'.join(ast.get_source_segment(source.read_text(encoding='utf-8'),n) for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in keep)
(HERE/'reference_allocation.py').write_text(code,encoding='utf-8')
(HERE/'inputs/PROVENANCE.json').write_text(json.dumps({'status':'PASS','source_npz_sha256':source_hashes,'allocation_source_sha256':sha(source),'original_count_tables_matched':200,'within_stratum_group_and_label_alignment':True,'group_counts_sha256':sha(HERE/'inputs/paired_group_counts.npz')},indent=2))
print('Prepared and reconciled 100 matched pairs, grouped by 10 cohort/batch strata.')
