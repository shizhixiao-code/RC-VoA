"""Run frozen-data reproductions on isolated copies; never overwrite archive results."""
from pathlib import Path
from datetime import datetime, timezone
import argparse, hashlib, json, os, platform, shutil, subprocess, sys, time, uuid, zipfile

ROOT = Path(__file__).resolve().parent

def verify_manifest():
    manifest = json.loads((ROOT / 'MANIFEST.json').read_text(encoding='utf8'))
    for name, expected in manifest['files'].items():
        p = ROOT / name
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest() != expected:
            raise RuntimeError(f'Release checksum mismatch: {name}')
    return len(manifest['files'])

def extract_safe(archive, destination):
    with zipfile.ZipFile(archive) as z:
        for item in z.infolist():
            p = (destination / item.filename).resolve()
            if not p.is_relative_to(destination.resolve()):
                raise ValueError(f'Unsafe archive member: {item.filename}')
        z.extractall(destination)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite', choices=['quick','core','quota','endpoints','verification-checks',
                        'verification-full','source-inclusion','bootstrap','standard'], default='quick')
    parser.add_argument('--work-dir', type=Path, default=ROOT/'runs')
    args = parser.parse_args()
    verified = verify_manifest()
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:6]
    work = args.work_dir.resolve()/f'{args.suite}_{stamp}'
    work.mkdir(parents=True, exist_ok=False)
    (work/'logs').mkdir()
    receipt = {'suite':args.suite,'started_utc':stamp,'status':'RUNNING','release_files_verified':verified,
               'python':platform.python_version(),'platform':platform.platform(),'steps':[]}
    for package in ['numpy','pandas','scipy']:
        from importlib.metadata import version
        receipt[package] = version(package)
    env = dict(os.environ, PYTHONIOENCODING='utf-8', OPENBLAS_NUM_THREADS='1')
    def save():
        (work/'REPLAY.json').write_text(json.dumps(receipt,indent=2),encoding='utf8')
    copied = {}
    def module(name):
        if name not in copied:
            dst=work/name
            shutil.copytree(ROOT/name,dst,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
            copied[name]=dst
            # Original workstation preservation is optional provenance, not a numerical input.
            if name=='endpoint_sensitivity':
                p=dst/'validation/PRESERVED_CURRENT_VERSION.json'
                data=json.loads(p.read_text())
                data['original_root']=str(work/'original_workspace_not_distributed')
                p.write_text(json.dumps(data,indent=2),encoding='utf8')
        return copied[name]
    def run(directory, script, *arguments):
        relative = str(directory.relative_to(work)) if directory.is_relative_to(work) else 'release/'+str(directory.relative_to(ROOT))
        step={'cwd':relative,'command':[script,*arguments]}
        receipt['steps'].append(step);save()
        print(f'RUN {step["cwd"]}: {script} {" ".join(arguments)}',flush=True)
        start=time.perf_counter()
        log=work/'logs'/f'{len(receipt["steps"]):02d}_{Path(script).stem}.txt'
        with log.open('w',encoding='utf8') as stream:
            completed=subprocess.run([sys.executable,str(directory/script),*arguments],cwd=directory,
                                     env=env,stdout=stream,stderr=subprocess.STDOUT)
        step.update(exit_code=completed.returncode,seconds=round(time.perf_counter()-start,3),log=str(log.relative_to(work)))
        save()
        if completed.returncode:
            print(log.read_text(encoding='utf8')[-6000:],flush=True)
            raise RuntimeError(f'Failed {script}; see {log}')
    def quick():
        c=module('core');s=module('source_data')
        run(c,'reproduce.py','--example')
        run(c,'replay_verification.py')
        run(s,'check_source_data.py')
        run(c/'label_sensitivity','check_table.py')
        run(c/'threshold_transport','check_transport_data.py')
        run(ROOT/'scripts','summarize_endpoints.py','--data-dir',str(ROOT/'endpoint_sensitivity'),
            '--output-dir',str(work/'endpoint_tables'))
    def core():
        c=module('core');run(c,'reproduce.py','--all','--check')
        run(c/'target_sensitivity','analyze_targets.py')
        run(c,'check_target_data.py')
    def quota():
        q=module('quota');run(q,'run.py','--replay');run(q,'check_prefix.py')
    def endpoints():
        e=module('endpoint_sensitivity');run(e/'code','run_analysis.py');run(e/'code','check_changed_case.py')
        run(ROOT/'scripts','summarize_endpoints.py','--data-dir',str(e),'--output-dir',str(work/'endpoint_tables'))
    def verification():
        e=module('core')/'exact_completion'
        run(e,'validate.py');run(e,'audit_boundaries_portable.py')
    def source_inclusion():
        dst=work/'expanded_source';dst.mkdir()
        extract_safe(ROOT/'archives/RCVoA_P1_Full_Results.zip',dst)
        p=dst/'rcvoa_p1_full_20261109/validation/PRESERVED_CURRENT_VERSION.json'
        data=json.loads(p.read_text());data['original_root']=str(work/'original_workspace_not_distributed')
        p.write_text(json.dumps(data,indent=2))
        code=dst/'rcvoa_p1_full_20261109/code'
        for mode in ['pipeline','frozen']:
            run(code,'run_full_analysis.py','--reference-mode',mode)
            run(code,'bootstrap_repair.py','--reference-mode',mode)
        run(code,'check_sccis_equivalence.py')
    try:
        if args.suite=='standard':
            for fn in [quick,core,quota,endpoints,verification]:fn()
        elif args.suite=='quick':quick()
        elif args.suite=='core':core()
        elif args.suite=='quota':quota()
        elif args.suite=='endpoints':endpoints()
        elif args.suite=='verification-checks':verification()
        elif args.suite=='verification-full':
            e=module('core')/'exact_completion'
            run(e,'run_experiment.py','full')
            run(e,'run_experiment.py','online','--orders','10')
            run(e,'summarize_portable.py')
        elif args.suite=='source-inclusion':source_inclusion()
        elif args.suite=='bootstrap':run(module('core')/'repair_uncertainty','run_bootstrap.py')
        receipt['status']='PASS'
    except Exception as error:
        receipt.update(status='FAIL',error=str(error));save();raise
    save();print(f'PASS: {work / "REPLAY.json"}',flush=True)

if __name__=='__main__':main()
