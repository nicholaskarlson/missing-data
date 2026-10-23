"""Raw-history and detailed risk-set audit, without analyzer imports."""
import csv,math

def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def close(a,b):return math.isfinite(float(a)) and math.isclose(float(a),float(b),abs_tol=1e-11,rel_tol=1e-11)
def compare_outputs(reference,actual,risk_reference,risk_actual,root):
 obs=sorted(read(root/'data/anchor_observed_registry.csv'),key=lambda r:r['student_id']);a,b=read(reference),read(actual)
 if len(b)!=len(obs) or len(a)!=len(obs) or [r['student_id'] for r in b]!=[r['student_id'] for r in obs] or a!=b:raise ValueError('Unit10 person identity/raw row mismatch')
 for row,raw in zip(b,obs):
  for k in raw:
   if k!='cohort' and row[k]!=raw[k]:raise ValueError('raw history changed '+k)
  group=('offered' if raw['offer']=='1' else 'control') if raw['responded_baseline']=='1' else 'nonrespondents'
  if row['group']!=group:raise ValueError('observation group changed')
  from decimal import Decimal
  if int(Decimal(raw['exit_month'])*1000)!=int(row['ticks']):raise ValueError('recorded clock changed')
 a,b=read(risk_reference),read(risk_actual)
 if len(a)!=len(b):raise ValueError('risk-set inventory')
 numeric=set(a[0])-{'group'}
 for x,y in zip(a,b):
  if set(x)!=set(y) or x['group']!=y['group'] or any(not close(x[k],y[k]) for k in numeric):raise ValueError('risk-set route disagreement')
 for g in ['pooled','control','offered','nonrespondents']:
  persons=obs if g=='pooled' else [r for r in read(actual) if r['group']==g]
  histories=[(int(Decimal(r['exit_month'])*1000),r['exit_type']) for r in persons];rr=[r for r in b if r['group']==g]
  times=sorted({t for t,e in histories})
  if [int(r['ticks']) for r in rr]!=times:raise ValueError('risk clock inventory')
  S=1.;W=T=0.;K=1.
  for r,t in zip(rr,times):
   n=sum(s>=t for s,e in histories);counts={e:sum(s==t and v==e for s,v in histories) for e in ['withdrawal','transfer','consent_withdrawn','none']}
   expected={'risk':n,'withdrawals':counts['withdrawal'],'transfers':counts['transfer'],'consent':counts['consent_withdrawn'],'administrative':counts['none'],'s_before':S,'increment_w':S*counts['withdrawal']/n,'increment_t':S*counts['transfer']/n}
   W+=expected['increment_w'];T+=expected['increment_t'];S*=1-(counts['withdrawal']+counts['transfer'])/n;K*=1-counts['withdrawal']/n;expected.update(withdrawal=W,transfer=T,event_free=S,withdrawal_only=1-K)
   if any(not close(r[k],v) for k,v in expected.items()) or not close(float(r['withdrawal'])+float(r['transfer'])+float(r['event_free']),1):raise ValueError('risk-set raw reconstruction')
 return {'status':'success','person_rows':len(obs),'riskset_rows':len(b),'unknown_outcomes_preserved':sum(r['outcome_known']=='0' and r['withdrawal_label']=='' for r in obs),'raw_risk_reconstruction':True}
