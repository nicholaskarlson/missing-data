"""Standalone companion checks and bounded, verified Unit 1 resume records."""
import argparse,csv,hashlib,json,math,os,subprocess,sys,time
from pathlib import Path,PurePosixPath
from compare_unit08 import compare_people as compare_weights
from compare_unit07 import compare_people as compare_likelihood
from compare_unit06 import compare_outputs as compare_joint
from compare_foundations import compare_people as compare_foundation_rows
from compare_unit01 import compare_people
from compare_unit09 import compare_people as compare_attrition
from compare_unit05 import compare_people as compare_completed
ROOT=Path(__file__).resolve().parents[1];M=ROOT/'open-materials/missing-data-reasoning-with-statistics'
CASES={'inventory':'observation.csv','mi':'mi_fixture.csv','withdrawal':'withdrawal_audit.csv','unit01':'unit01.csv','regmi':'regmi.csv','unit05':'unit05.csv','unit09':'unit09.csv','foundations':'foundations.csv','unit06':'unit06.csv','unit07':'unit07.csv','unit08':'unit08.csv','unit10':'unit10.csv','unit11':'unit11.csv','unit12':'unit12.csv','extension_a':'extension_a.csv','extension_b':'extension_b.csv'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def safe(n):
 if not n or PurePosixPath(n).is_absolute() or '..' in PurePosixPath(n).parts or '\\' in n:raise ValueError('unsafe export path')
 p=ROOT/n
 if not p.resolve().is_relative_to(ROOT) or any(x.is_symlink() for x in [p,*p.parents]):raise ValueError('symlink export path')
 return p

def boundary():
 pairs=[line.split('  ',1) for line in (ROOT/'SHA256SUMS').read_text().splitlines()];expected={n:h for h,n in pairs}
 if len(expected)!=len(pairs):raise ValueError('duplicate checksum path')
 actual={str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and not any(x in {'.git','build','__pycache__','.venv'} for x in p.relative_to(ROOT).parts)}
 if actual!=set(expected)|{'SHA256SUMS'}:raise ValueError('export inventory differs')
 manifest=json.loads((ROOT/'EXPORT_MANIFEST.json').read_text())
 if set(manifest['files'])|{'SHA256SUMS','EXPORT_MANIFEST.json'}!=actual:raise ValueError('positive allowlist differs')
 for n,h in expected.items():
  if sha(safe(n))!=h:raise ValueError('changed export file '+n)
  if any(x in {'manuscript','templates','hidden','author','audit','reveal'} for x in PurePosixPath(n).parts):raise ValueError('private category in ordinary export')
 for n,h in manifest['files'].items():
  if sha(safe(n))!=h:raise ValueError('source identity differs '+n)
 catalog=json.loads((ROOT/'site/catalog.json').read_text())
 if catalog['implemented_units']!=12 or [u['number'] for u in catalog['units']]!=list(range(1,13)) or any(u['status']!='implemented' for u in catalog['units']):raise ValueError('unit availability misrepresented')
 if catalog.get('implemented_numerical_cases')!=list(CASES) or catalog.get('numerical_case_count')!=len(CASES):raise ValueError('numerical case catalog differs from executable verifier')
 if catalog.get('repository')!='https://github.com/nicholaskarlson/missing-data' or catalog.get('release')!='v0.1.0' or manifest.get('release')!=catalog['release'] or manifest.get('repository')!='nicholaskarlson/missing-data':raise ValueError('release identity differs')
 from html.parser import HTMLParser
 class Links(HTMLParser):
  def __init__(self):super().__init__();self.links=[]
  def handle_starttag(self,tag,attrs):
   for k,v in attrs:
    if k in {'href','src'} and v:self.links.append(v)
 for page in (ROOT/'site').glob('*.html'):
  parser=Links();parser.feed(page.read_text())
  for link in parser.links:
   if link.startswith(('#','https://','http://','mailto:')):continue
   target=(page.parent/link.split('#',1)[0]).resolve()
   if not target.is_relative_to(ROOT) or not target.is_file():raise ValueError('broken site link '+link)
 if catalog['units'][0].get('page')!='unit01.html' or not (ROOT/'site/unit01.html').is_file():raise ValueError('Unit 1 guide missing')
 if catalog.get('method_pilots')!=[{'unit':5,'case_id':'B8-U05-CAL-001','status':'calibration_only','page':'regmi.html','command':'python scripts/check_regmi.py --base-r --outdir build/regmi-pilot'}] or not (ROOT/'site/regmi.html').is_file():raise ValueError('calibration scope/guide missing')
 if catalog['units'][4].get('page')!='unit05.html' or catalog['units'][4].get('model_audit')!='conditional_independence_warning' or 'Model audit: warning.' not in (ROOT/'site/unit05.html').read_text():raise ValueError('Unit 5 guide/warning missing')
 for n in [2,3,4]:
  if catalog['units'][n-1].get('page')!='unit%02d.html'%n or not (ROOT/('site/unit%02d.html'%n)).is_file():raise ValueError('foundation guide missing')
 if catalog['units'][8].get('page')!='unit09.html' or catalog['units'][8].get('inference')!='support_and_conservative_randomization' or 'randomization' not in (ROOT/'site/unit09.html').read_text():raise ValueError('Unit 9 guide/inference missing')
 if catalog['units'][5].get('page')!='unit06.html' or catalog['units'][5].get('target')!='B8-D06' or not (ROOT/'site/unit06.html').is_file():raise ValueError('Unit6 diagnostic scope missing')
 if catalog['units'][6].get('page')!='unit07.html' or catalog['units'][6].get('target')!='B8-M07' or not (ROOT/'site/unit07.html').is_file():raise ValueError('Unit7 model probability scope missing')
 if catalog['units'][7].get('page')!='unit08.html' or catalog['units'][7].get('target')!='B8-T01' or not (ROOT/'site/unit08.html').is_file():raise ValueError('Unit8 weighting scope missing')
 if catalog['units'][10].get('page')!='unit11.html' or catalog['units'][10].get('inference')!='mean_sensitivity_without_interval' or not (ROOT/'site/unit11.html').is_file():raise ValueError('Unit11 sensitivity guide/scope missing')
 if catalog['units'][11].get('page')!='unit12.html' or catalog['units'][11].get('target')!='B8-D12' or catalog['units'][11].get('inference')!='conditional_form_design_conservative_interval' or not (ROOT/'site/unit12.html').is_file():raise ValueError('Unit12 conditional form scope missing')
 for extension in catalog.get('extensions',[]):
  if extension.get('case') not in {'extension_a','extension_b'} or extension.get('status')!='implemented' or not (ROOT/'site'/extension['page']).is_file():raise ValueError('extension guide missing')
 if len(catalog.get('extensions',[]))!=2 or {e.get('id') for e in catalog.get('extensions',[])}!={'A','B'}:raise ValueError('two extension guides required')
 if not manifest.get('fresh_history_required') or manifest.get('reveal_included') is not False:raise ValueError('ordinary publication boundary differs')
 print('MISSING_BOUNDARY_OK files='+str(len(actual))+' implemented_units=12',flush=True)
 return manifest

def values(p):
 with p.open() as f:rows=list(csv.DictReader(f))
 if len({r['key'] for r in rows})!=len(rows):raise ValueError('duplicate output key')
 return {r['key']:float(r['value']) for r in rows}

def execute(out,cases,base_r=False):
 boundary();import numpy as np
 if np.__version__!='1.26.4':raise ValueError('NumPy 1.26.4 required')
 out=out.resolve();out.mkdir(parents=True,exist_ok=False);routes=['numpy','stdlib']+(['base_r'] if base_r else []);report={}
 for case in cases:
  reference=values(M/'expected-results'/CASES[case]);report[case]={}
  for route in routes:
   output=out/(case+'_'+route+'.csv');script=M/'source'/(case+'_'+route+('.R' if route=='base_r' else '.py'));command=(['Rscript','--vanilla'] if route=='base_r' else [sys.executable])+[str(script),str(M),str(output)]
   r=subprocess.run(command,capture_output=True,text=True,timeout=45);(out/(case+'_'+route+'.log')).write_text(r.stdout+r.stderr)
   if r.returncode or case.upper()+'_'+route.upper()+'_OK' not in r.stdout:raise ValueError('route failed '+case+'/'+route+'\n'+r.stdout+r.stderr)
   actual=values(output)
   if set(actual)!=set(reference):raise ValueError('output key inventory differs')
   for k,v in actual.items():
    if not math.isfinite(v) or not math.isclose(v,reference[k],abs_tol=1e-11,rel_tol=1e-11):raise ValueError('numerical disagreement '+case+'/'+k)
    if (case=='inventory' or k.startswith('pattern/')) and v!=int(v):raise ValueError('noninteger count')
   if case=='regmi':
    with (M/'data/regmi_calibration.csv').open() as stream:original={r['row_id']:float(r['y']) for r in csv.DictReader(stream) if r['y']!=''}
    for j in range(1,21):
     for sid,value in original.items():
      if actual[f'imputation/{j}/completed/{sid}']!=value:raise ValueError('observed calibration value changed')
   record={'status':'success','output_sha256':sha(output)}
   if case in {'extension_a','extension_b'}:
    from compare_extensions import compare_outputs as compare_extension
    record['rows']=compare_extension(case,M/'expected-results'/(case+'.csv'),output,M)
   if case=='mi':
    with (M/'data/mi_fixture.csv').open() as stream:original={r['row_id']:float(r['y']) for r in csv.DictReader(stream) if r['y']!=''}
    if any(actual[f'imputation/{j}/completed/{sid}']!=v for j in range(1,9) for sid,v in original.items()):raise ValueError('observed fixture value changed')
   if case=='unit01':record['people']=compare_people(M/'expected-results/unit01_people.csv',out/(case+'_'+route+'_people.csv'))
   if case=='unit05':record['people']=compare_completed(M/'expected-results/unit05_people.csv',out/(case+'_'+route+'_people.csv'),M)
   if case=='unit09':record['people']=compare_attrition(M/'expected-results/unit09_people.csv',out/(case+'_'+route+'_people.csv'),M)
   if case=='foundations':record['people']=compare_foundation_rows(M/'expected-results/foundations_people.csv',out/(case+'_'+route+'_people.csv'),M)
   if case=='unit08':
    record['rows']=compare_weights(M/'expected-results/unit08_people.csv',out/(case+'_'+route+'_people.csv'),M)
    if not math.isclose(record['rows']['mean_reconstructed_from_rows'],actual['estimate/poststratified_mean'],abs_tol=1e-11,rel_tol=1e-11):raise ValueError('Unit8 row estimate differs')
   if case=='unit12':
    from compare_unit12 import compare_outputs as compare_forms
    record['rows']=compare_forms(M/'expected-results/unit12.csv',out/(case+'_'+route+'.csv'),M)
   if case=='unit11':
    from compare_unit11 import compare_outputs as compare_sensitivity
    record['rows']=compare_sensitivity(M/'expected-results/unit11.csv',out/(case+'_'+route+'.csv'),M)
   if case=='unit10':
    from compare_unit10 import compare_outputs as compare_events
    record['rows']=compare_events(M/'expected-results/unit10_people.csv',out/(case+'_'+route+'_people.csv'),M/'expected-results/unit10_risksets.csv',out/(case+'_'+route+'_risksets.csv'),M)
   if case=='unit07':
    record['rows']=compare_likelihood(M/'expected-results/unit07_people.csv',out/(case+'_'+route+'_people.csv'),M)
    if not math.isclose(record['rows']['plugin_prediction'],actual['finite/plugin_prediction'],abs_tol=1e-11,rel_tol=1e-11):raise ValueError('Unit7 person estimate disagreement')
   if case=='unit06':
    record['rows']=compare_joint(M/'expected-results/unit06_people.csv',out/(case+'_'+route+'_people.csv'),M/'expected-results/unit06_completed.csv',out/(case+'_'+route+'_completed.csv'),M)
    for d,n in record['rows']['completion_event_totals'].items():
     if actual['draw/'+d+'/Q']!=n/1600:raise ValueError('completed estimate mismatch')
   report[case][route]=record
 (out/'report.json').write_text(json.dumps({'status':'success','cases':report,'implemented_units':12,'unit05_model_audit':'warning' if 'unit05' in cases else 'not_run','numpy':np.__version__,'python':sys.version,'remote_ci_claimed':False},indent=2)+'\n')
 print('MISSING_NUMERICAL_OK cases='+str(len(cases))+' routes='+str(len(routes))+' implemented_units=12',flush=True)
 return report


def valid_resume(rec,fingerprint,base_r):
 try:
  if rec.get('status')!='success' or rec.get('fingerprint')!=fingerprint:return False
  routes=['numpy','stdlib']+(['base_r'] if base_r else [])
  if rec['routes']!=routes:return False
  attempt=safe(rec['output_folder']);needed={'summary.csv','people.csv','report.json'}|{name for r in routes for name in ['unit01_'+r+'.csv','unit01_'+r+'_people.csv','unit01_'+r+'.log']}
  if {p.name for p in attempt.iterdir() if p.is_file()}!=needed:return False
  if set(rec['output_sha256'])!={str(attempt.relative_to(ROOT)/n) for n in needed}:return False
  if not all(sha(safe(n))==h for n,h in rec['output_sha256'].items()):return False
  reference=values(M/'expected-results/unit01.csv')
  for route in routes:
   actual=values(attempt/('unit01_'+route+'.csv'))
   if set(actual)!=set(reference) or any(not math.isfinite(v) or not math.isclose(v,reference[k],abs_tol=1e-11,rel_tol=1e-11) for k,v in actual.items()):return False
   if 'UNIT01_'+route.upper()+'_OK' not in (attempt/('unit01_'+route+'.log')).read_text():return False
   compare_people(M/'expected-results/unit01_people.csv',attempt/('unit01_'+route+'_people.csv'))
  if (attempt/'summary.csv').read_bytes()!=(attempt/'unit01_numpy.csv').read_bytes() or (attempt/'people.csv').read_bytes()!=(attempt/'unit01_numpy_people.csv').read_bytes():return False
  return True
 except (OSError,ValueError,KeyError,TypeError):return False

def unit_cli():
 p=argparse.ArgumentParser();p.add_argument('--numeric',action='store_true');p.add_argument('--base-r',action='store_true');a=p.parse_args();boundary()
 if not(a.numeric or a.base_r):return
 import fcntl,numpy as np
 folder=ROOT/'build/unit01-audit';folder.mkdir(parents=True,exist_ok=True)
 with (folder/'runner.lock').open('a+') as lock:
  try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  except BlockingIOError:raise ValueError('another Unit 1 writer is active')
  version={'python':sys.version,'numpy':np.__version__}
  if a.base_r:
   r=subprocess.run(['Rscript','--vanilla','-e','cat(as.character(getRversion()))'],capture_output=True,text=True,check=True,timeout=20);version['R']=r.stdout
  manifest=json.loads((ROOT/'EXPORT_MANIFEST.json').read_text());names=[n for n in manifest['files'] if any(x in n for x in ['source/unit01_','data/anchor_participants.csv','data/anchor_baseline.csv','model-specs/unit01','expected-results/unit01','scripts/compare_unit01.py','scripts/verify_companion.py','scripts/check_unit01.py'])]
  fingerprint=hashlib.sha256(json.dumps({'files':{n:sha(safe(n)) for n in names},'runtime':version,'base_r':a.base_r},sort_keys=True).encode()).hexdigest();current=folder/'current.json'
  if current.exists():
   try:rec=json.loads(current.read_text())
   except (OSError,ValueError):rec={}
   if valid_resume(rec,fingerprint,a.base_r):
    print('MISSING_UNIT01_RESUME sources_runtime_outputs_verified=true output_folder='+rec['output_folder']);return
  attempt=folder/('attempt-'+str(time.time_ns()));execute(attempt,['unit01'],a.base_r)
  shutil=__import__('shutil');shutil.copyfile(attempt/'unit01_numpy.csv',attempt/'summary.csv');shutil.copyfile(attempt/'unit01_numpy_people.csv',attempt/'people.csv')
  record={'status':'success','case':'unit01','fingerprint':fingerprint,'versions':version,'routes':['numpy','stdlib']+(['base_r'] if a.base_r else []),'output_folder':str(attempt.relative_to(ROOT)),'output_sha256':{str(f.relative_to(ROOT)):sha(f) for f in attempt.iterdir() if f.is_file()},'remote_ci_claimed':False}
  tmp=current.with_suffix('.tmp');tmp.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n');tmp.replace(current)
  v=values(attempt/'summary.csv');print('MISSING_UNIT01_OK complete_people='+str(int(v['complete_people']))+' complete_mean='+format(v['complete_mean'],'.3f')+' item_lower='+format(v['item_lower'],'.5f')+' item_upper='+format(v['item_upper'],'.5f')+' output_folder='+record['output_folder'])

def main():
 p=argparse.ArgumentParser();p.add_argument('--base-r',action='store_true');p.add_argument('--boundary-only',action='store_true');p.add_argument('--case',choices=['all',*CASES],default='all');p.add_argument('--outdir',type=Path,default=ROOT/'build/check');a=p.parse_args();boundary()
 if not a.boundary_only:execute(a.outdir,list(CASES) if a.case=='all' else [a.case],a.base_r)
if __name__=='__main__':main()
