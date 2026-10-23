"""Copy only verified ordinary export files into a new Pages artifact tree."""
import argparse,hashlib,json,shutil
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit,unquote
from verify_companion import ROOT,boundary,safe
def links(root):
 class Scan(HTMLParser):
  def __init__(self):super().__init__();self.refs=[];self.ids=set()
  def handle_starttag(self,tag,attrs):
   for k,v in attrs:
    if k in {'href','src'} and v:self.refs.append(v)
    if k=='id' and v:self.ids.add(v)
 scans={}
 for p in root.rglob('*.html'):
  s=Scan();s.feed(p.read_text());scans[p.resolve()]=s
 for page,s in scans.items():
  for ref in s.refs:
   u=urlsplit(ref)
   if u.scheme or u.netloc:continue
   # Root-relative URLs would escape the GitHub project prefix.
   if u.path.startswith('/'):raise ValueError('project-prefix escaping link '+ref)
   target=(page.parent/unquote(u.path)).resolve() if u.path else page
   if not target.is_relative_to(root.resolve()) or not target.is_file():raise ValueError('broken deployed link '+str(page)+' '+ref)
   if u.fragment and target.suffix=='.html' and unquote(u.fragment) not in scans[target].ids:raise ValueError('missing deployed anchor '+ref)
 return len(scans)
def build(out):
 manifest=boundary();out=out.expanduser().absolute()
 if out.exists() or any(p.is_symlink() for p in [out,*out.parents]):raise ValueError('new nonsymlink site directory required')
 out.mkdir(parents=True)
 selected={n:h for n,h in manifest['files'].items() if not Path(n).parts[0].startswith('.')}
 for n in selected:
  dest=out/n;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(safe(n),dest)
 pages=links(out)
 (out/'SITE_MANIFEST.json').write_text(json.dumps({'kind':'Pages projection of verified numerical export','files':selected,'release':manifest['release'],'source_export_manifest_sha256':hashlib.sha256((ROOT/'EXPORT_MANIFEST.json').read_bytes()).hexdigest(),'private_material_included':False},indent=2,sort_keys=True)+'\n')
 print('MISSING_SITE_OK pages='+str(pages)+' files='+str(len(selected))+' project_prefix=/missing-data/',flush=True)
 return selected
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);build(p.parse_args().outdir)
