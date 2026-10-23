"""Array lattice integration of the declared profile joint model."""
import csv,json,sys
from pathlib import Path
import numpy as np

def calculate(root):
 root=Path(root);s=json.loads((root/'model-specs/extension_b.json').read_text());items=s['items']
 if s['n']!=6 or s['prior_per_profile']!=1 or s['reverse_item']!='st4' or items!=['st'+str(i) for i in range(1,6)]:raise ValueError('B contract')
 def read(n):
  with (root/n).open() as f:return list(csv.DictReader(f))
 rows=read('data/extension_b.csv');p=read('model-specs/extension_b_profiles.csv')
 if [r['person'] for r in rows]!=['B'+str(i) for i in range(1,7)] or [r['profile'] for r in p]!=list('ABCD'):raise ValueError('B keys')
 if any(v[j] not in {'','1','2','3','4','5'} for v in rows for j in items) or any(v[j] not in {'1','2','3','4','5'} for v in p for j in items):raise ValueError('B values')
 obs=np.array([[float(v[j]) if v[j] else np.nan for j in items] for v in rows]);raw=np.array([[int(v[j]) for j in items] for v in p]);keyed=raw.copy();
 if len(np.unique(raw,axis=0))!=4:raise ValueError('B duplicate profiles')
 keyed[:,3]=6-keyed[:,3];score=keyed.mean(axis=1)
 allowed=np.all(np.isnan(obs[:,None,:])|(obs[:,None,:]==raw[None,:,:]),axis=2)
 if not np.all(allowed.any(axis=1)):raise ValueError('B impossible profile')
 counts=np.sum(~np.isnan(obs),axis=1);known=np.flatnonzero(counts==5);unk=np.flatnonzero(counts<5);alpha=np.ones(4)
 for i in known:alpha[np.argmax(allowed[i])]+=1
 lattice=np.indices((4,)*len(unk)).reshape(len(unk),-1).T if len(unk) else np.empty((1,0),int)
 z=np.zeros((len(lattice),6),int);z[:,unk]=lattice
 for i in known:z[:,i]=np.argmax(allowed[i])
 cs=np.stack([(lattice==k).sum(axis=1) for k in range(4)],axis=1);w=np.ones(len(lattice))
 for k in range(4):
  for j in range(len(unk)):w*=np.where(cs[:,k]>j,alpha[k]+j,1)
 w/=np.prod(alpha.sum()+np.arange(len(unk)));Q=score[z].mean(axis=1);keep=np.all(allowed[np.arange(6)[None,:],z],axis=1)
 evidence=float(w[keep].sum());pi=np.where(keep,w,0)/evidence;ps=w/w.sum()
 out={'data/n':6,'data/full_scores':len(known),'data/partial_people':np.sum((counts>0)&(counts<5)),'data/observed_item_cells':counts.sum(),'data/discarded_partial_cells':counts[unk].sum(),'data/state_count':len(lattice)}
 for route,prob,e in [('item',pi,evidence),('score',ps,w.sum())]:
  mu=prob@Q;out.update({route+'/evidence':e,route+'/mean':mu,route+'/variance':prob@((Q-mu)**2)})
 val=obs.copy();val[:,3]=6-val[:,3];sums=np.nansum(val,axis=1);lo=(sums+5-counts)/5;hi=(sums+5*(5-counts))/5;out.update(logical_lower=lo.mean(),logical_upper=hi.mean(),available_mean_n=np.sum(counts>0),available_mean=np.mean(sums[counts>0]/counts[counts>0]))
 people=[]
 for i,v in enumerate(rows):
  row=dict(person=v['person'],observed_items=counts[i],logical_lower=lo[i],logical_upper=hi[i])
  for route,prob in [('item',pi),('score',ps)]:mu=prob@score[z[:,i]];row[route+'_mean']=mu;row[route+'_variance']=prob@((score[z[:,i]]-mu)**2)
  people.append(row)
 states=[dict(assignment=''.join('ABCD'[j] for j in l),Q=q,item_probability=a,score_probability=b) for l,q,a,b in zip(lattice,Q,pi,ps)]
 return out,people,states

def write(path,rows,fields):
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
  for row in rows:w.writerow({k:v if isinstance(v,str) else format(float(v),'.17g') for k,v in row.items()})
def main():
 root,out=map(Path,sys.argv[1:]);v,p,t=calculate(root);write(out,[{'key':k,'value':v[k]} for k in sorted(v)],['key','value']);write(out.with_name(out.stem+'_people.csv'),p,['person','observed_items','logical_lower','logical_upper','item_mean','item_variance','score_mean','score_variance']);write(out.with_name(out.stem+'_states.csv'),t,['assignment','Q','item_probability','score_probability']);print('EXTENSION_B_NUMPY_OK')
if __name__=='__main__':main()
