"""Reconstruct rows, missingness ledger and every teaching allocation from raw inputs."""
import csv,math,itertools
from fractions import Fraction as F
ITEMS=['wb1','wb2','wb3','wb4','wb5'];AUDIT=ITEMS+['st1','st2','st3','st4','st5']
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def close(a,b):return math.isfinite(float(a)) and math.isclose(float(a),float(b),abs_tol=1e-11,rel_tol=1e-11)
def require(x,m):
 if not x:raise ValueError('Unit12 '+m)
def key(r,k):return None if r[k]=='' else 6-int(r[k]) if k in ['wb4','st4'] else int(r[k])
def compare_outputs(reference,actual,root):
 data={}
 for suffix in ['people','cells','masks']:
  a=read(reference.with_name(reference.stem+'_'+suffix+'.csv'));b=read(actual.with_name(actual.stem+'_'+suffix+'.csv'));require(len(a)==len(b) and set(a[0])==set(b[0]),'detail inventory')
  exact={'student_id','college_id','offer','responded6','form','item','fixture','forms','full_score',*AUDIT}
  for x,y in zip(a,b):
   for k in x:require(x[k]==y[k] if k in exact or x[k]=='' else close(x[k],y[k]),'reference detail '+k)
  data[suffix]=b
 p=data['people'];c=data['cells'];masks=data['masks'];base=sorted(read(root/'data/anchor_baseline.csv'),key=lambda r:r['student_id']);six={r['student_id']:r for r in read(root/'data/anchor_followup.csv') if r['wave_month']=='6'}
 require(len(p)==1219 and len(c)==40 and len(masks)==243,'fixed detailed counts');hlo=[];hhi=[];full=[[],[]]
 for b,row in zip(base,p):
  r=six[b['student_id']];require(all(row[k]==b[k] for k in ['student_id','college_id','offer']) and row['responded6']==r['responded'] and row['form']==r['form'] and all(row[k]==r[k] for k in AUDIT),'raw answers or blanks changed');yes=r['responded']=='1';planned={f'{block}{int(r["form"])+2}' for block in ['wb','st']} if yes else set();x={k:key(r,k) for k in AUDIT}
  for block in ['wb','st']:
   items=[k for k in AUDIT if k.startswith(block)];expect={'observed':sum(x[k] is not None for k in items),'planned':int(yes),'wave':0 if yes else 5,'incidental':sum(x[k] is None and k not in planned for k in items) if yes else 0}
   require(all(close(row[block+'_'+k],v) for k,v in expect.items()),'person partition')
  ks=sum(x[k] for k in ITEMS if x[k] is not None);mi=sum(x[k] is None for k in ITEMS);bounds=[F(ks+mi,5),F(ks+5*mi,5)];require(row['full_score']=='' and close(row['keyed_wb_sum'],ks) and close(row['full_lower'],bounds[0]) and close(row['full_upper'],bounds[1]),'logical score support')
  for i,v in enumerate(bounds):full[i].append(v)
  if yes:
   lo=F(0);hi=F(0);known=F(0)
   for j,k in enumerate(ITEMS):
    if k in planned:continue
    w=F(1) if j<2 else F(3,2)
    if x[k] is None:lo+=w;hi+=5*w
    else:lo+=w*x[k];hi+=w*x[k];known+=w*x[k]
   require(close(row['H_lower'],lo/5) and close(row['H_upper'],hi/5) and close(row['H_known_total'],known),'weighted contribution bounds');hlo.append(lo/5);hhi.append(hi/5)
  else:require(all(row[k]=='' for k in ['H_known_total','H_lower','H_upper']),'nonrespondent in conditional design')
 require([(r['form'],r['item']) for r in c]==[(str(f),k) for f in range(4) for k in AUDIT],'ledger keys')
 for cell in c:
  form,k=cell['form'],cell['item'];rows=[r for r in p if r['form']==form];observed=[key(r,k) for r in rows if r[k]];wave=len(rows) if form=='0' else 0;planned=len(rows) if form!='0' and k[-1]==str(int(form)+2) else 0;expected=dict(roster_n=len(rows),observed_n=len(observed),planned_n=planned,wave_n=wave,incidental_n=len(rows)-len(observed)-wave-planned,keyed_sum=sum(observed));require(all(close(cell[k],v) for k,v in expected.items()),'cell reconstruction')
 toy=read(root/'model-specs/unit12_teaching.csv');wanted=[];toy_variance={}
 for fixture in ['parallel','mixed','unequal_blocks']:
  rows=sorted([r for r in toy if r['fixture']==fixture],key=lambda r:r['id']);target=F(sum(key(r,k) for r in rows for k in ITEMS),20);hh=[];vv=[]
  for forms in itertools.product('123',repeat=4):
   h=F(0);nv=F(0);v=F(0)
   for r,f in zip(rows,forms):
    common=key(r,'wb1')+key(r,'wb2');pair=[key(r,k) for k in ITEMS[2:] if k!='wb'+str(int(f)+2)];h+=F(2*common+3*sum(pair),40);nv+=F(common+sum(pair),16);v+=F(3*(pair[0]-pair[1])**2,1600)
   wanted.append(dict(fixture=fixture,forms=''.join(forms),target=target,H=h,naive=nv,Vhat=v));hh.append(h);vv.append(v)
  mean=sum(hh)/81;var=sum((h-mean)**2 for h in hh)/81;require(mean==target and sum(vv)/81==var,'enumerated unbiasedness/variance');toy_variance[fixture]=float(var)
 for a,b in zip(wanted,masks):require(all(a[k]==b[k] if k in ['fixture','forms'] else close(a[k],b[k]) for k in a),'raw allocation reconstruction')
 s={r['key']:float(r['value']) for r in read(actual)};radius=math.sqrt(8/(25*956*.05));require(close(s['H/lower'],sum(hlo)/956) and close(s['H/upper'],sum(hhi)/956) and close(s['design/lower'],max(1,sum(hlo)/956-radius)) and close(s['design/upper'],min(5,sum(hhi)/956+radius)),'summary interval reconstruction')
 for pop,ix in [('randomized',list(range(1219))),('respondents',[i for i,r in enumerate(p) if r['responded6']=='1'])]:
  for j,side in enumerate(['lower','upper']):require(close(s['support/'+pop+'/'+side],sum(full[j][i] for i in ix)/len(ix)),'summary logical support')
 return dict(status='success',people_rows=1219,item_cells=40,mask_rows=243,raw_reconstruction=True,toy_variance=toy_variance)
