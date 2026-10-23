"""Exact rational integration over full ordered profile assignments."""
import csv,itertools,json,math,sys
from pathlib import Path
from fractions import Fraction as F

def calculate(root):
 root=Path(root);s=json.loads((root/'model-specs/extension_b.json').read_text());items=s['items']
 if s['n']!=6 or s['prior_per_profile']!=1 or s['reverse_item']!='st4' or items!=['st'+str(i) for i in range(1,6)]:raise ValueError('B contract')
 def read(n):
  with (root/n).open() as f:return list(csv.DictReader(f))
 r=read('data/extension_b.csv');pr=read('model-specs/extension_b_profiles.csv')
 if [v['person'] for v in r]!=['B'+str(i) for i in range(1,7)] or [v['profile'] for v in pr]!=list('ABCD'):raise ValueError('B keys')
 if any(v[j] not in {'','1','2','3','4','5'} for v in r for j in items) or any(v[j] not in {'1','2','3','4','5'} for v in pr for j in items):raise ValueError('B values')
 raw=[[int(v[j]) for j in items] for v in pr]
 if len(set(tuple(v) for v in raw))!=4:raise ValueError('B duplicate profiles')
 scores=[F(sum(x if j!=3 else 6-x for j,x in enumerate(row)),5) for row in raw]
 allowed=[[k for k,row in enumerate(raw) if all(v[j]=='' or int(v[j])==row[t] for t,j in enumerate(items))] for v in r]
 if any(not a for a in allowed):raise ValueError('B impossible profile')
 known=[i for i,v in enumerate(r) if all(v[j]!='' for j in items)];unknown=[i for i in range(6) if i not in known];alpha=[1]*4
 for i in known:alpha[allowed[i][0]]+=1
 states=[];weights={'item':[],'score':[]};totals=[];assignments=[]
 def rise(x,n):return math.prod(range(x,x+n))
 for assignment in itertools.product(range(4),repeat=len(unknown)):
  count=[assignment.count(k) for k in range(4)];w=F(math.prod(rise(alpha[k],count[k]) for k in range(4)),rise(sum(alpha),len(unknown)))
  z=[allowed[i][0] if i in known else assignment[unknown.index(i)] for i in range(6)];Q=sum(scores[k] for k in z)/6;ok=all(z[i] in allowed[i] for i in unknown)
  weights['item'].append(w if ok else F(0));weights['score'].append(w);totals.append(Q);assignments.append(z)
 out={'data/n':6,'data/full_scores':len(known),'data/partial_people':sum(0<sum(v[j]!='' for j in items)<5 for v in r),'data/observed_item_cells':sum(v[j]!='' for v in r for j in items),'data/discarded_partial_cells':sum(v[j]!='' for i,v in enumerate(r) if i not in known for j in items),'data/state_count':len(totals)}
 for route,ww in weights.items():
  evidence=sum(ww);prob=[w/evidence for w in ww];mu=sum(w*q for w,q in zip(prob,totals));out.update({route+'/evidence':evidence,route+'/mean':mu,route+'/variance':sum(w*(q-mu)**2 for w,q in zip(prob,totals))});weights[route]=prob
 people=[];lower=upper=F(0);avail=[]
 for i,v in enumerate(r):
  keyed=[int(v[j]) if t!=3 else 6-int(v[j]) for t,j in enumerate(items) if v[j]!=''];n=len(keyed);lo=F(sum(keyed)+5-n,5);hi=F(sum(keyed)+5*(5-n),5);lower+=lo/6;upper+=hi/6
  if n:avail.append(F(sum(keyed),n))
  row={'person':v['person'],'observed_items':n,'logical_lower':lo,'logical_upper':hi}
  for route,prob in weights.items():
   mu=sum(w*scores[z[i]] for w,z in zip(prob,assignments));row[route+'_mean']=mu;row[route+'_variance']=sum(w*(scores[z[i]]-mu)**2 for w,z in zip(prob,assignments))
  people.append(row)
 out.update(logical_lower=lower,logical_upper=upper,available_mean_n=len(avail),available_mean=sum(avail)/len(avail))
 for j,(Q,z) in enumerate(zip(totals,assignments)):states.append({'assignment':''.join('ABCD'[z[i]] for i in unknown),'Q':Q,'item_probability':weights['item'][j],'score_probability':weights['score'][j]})
 return out,people,states

def write(path,rows,fields):
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
  for row in rows:w.writerow({k:v if isinstance(v,str) else format(float(v),'.17g') for k,v in row.items()})
def main():
 root,out=map(Path,sys.argv[1:]);v,p,t=calculate(root);write(out,[{'key':k,'value':v[k]} for k in sorted(v)],['key','value']);write(out.with_name(out.stem+'_people.csv'),p,['person','observed_items','logical_lower','logical_upper','item_mean','item_variance','score_mean','score_variance']);write(out.with_name(out.stem+'_states.csv'),t,['assignment','Q','item_probability','score_probability']);print('EXTENSION_B_STDLIB_OK')
if __name__=='__main__':main()
