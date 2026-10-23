"""Always-run public status; unknown changes require numerical verification."""
import json,os,subprocess
from pathlib import Path
def presentation(n):
 return n in {'SHA256SUMS','CITATION.cff','index.html'} or n.endswith('.md') or (n.startswith('site/') and n.endswith(('.html','.css','.json')))
def projection(manifest):
 return {'files':{n:h for n,h in manifest['files'].items() if not presentation(n)},'identity':{k:manifest[k] for k in ['repository','release','completed_units','implemented_extensions','numerical_cases','fresh_history_required','reveal_included']}}
def scientific(names):
 for n in names:
  if presentation(n) or n=='EXPORT_MANIFEST.json':continue
  return True
 return False
def main():
 p=os.environ.get('GITHUB_EVENT_PATH');event=json.loads(Path(p).read_text()) if p else {};base=event.get('pull_request',{}).get('base',{}).get('sha') or event.get('before');head=os.environ.get('GITHUB_SHA','HEAD');required=True
 if base and set(base)!={'0'}:
  names=subprocess.check_output(['git','diff','--name-only',base,head],text=True,timeout=20).splitlines();required=scientific(names)
  if not required:
   try:
    before=json.loads(subprocess.check_output(['git','show',base+':EXPORT_MANIFEST.json'],text=True,timeout=20));after=json.loads(Path('EXPORT_MANIFEST.json').read_text());required=projection(before)!=projection(after)
   except (KeyError,ValueError,OSError,subprocess.SubprocessError):required=True
 value='true' if required else 'false';print('NUMERICAL_REQUIRED='+value)
 if os.environ.get('GITHUB_OUTPUT'):
  with open(os.environ['GITHUB_OUTPUT'],'a') as f:f.write('numerical='+value+'\n')
if __name__=='__main__':main()
