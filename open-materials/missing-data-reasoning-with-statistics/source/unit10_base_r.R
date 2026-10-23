# Independent base R three-state transition-matrix products.
args<-commandArgs(TRUE);root<-args[1];output<-args[2]
read<-function(n) read.csv(file.path(root,n),colClasses='character',na.strings=NULL,check.names=FALSE)
need<-function(x,msg) if(!isTRUE(x)) stop(msg)
ticks<-function(s) {
 need(length(s)==1 && grepl('^[0-9]+(\\.[0-9]{1,3})?$',s),'invalid recorded clock')
 a<-strsplit(s,'.',fixed=TRUE)[[1]];b<-if(length(a)==2) paste0(a[2],strrep('0',3-nchar(a[2]))) else '0';t<-as.integer(a[1])*1000L+as.integer(b)
 need(!is.na(t)&&t>=0&&t<=12000,'clock outside horizon');t
}
s<-read('model-specs/unit10_r.csv');need(!anyDuplicated(s$name),'spec duplicate')
spec<-setNames(s$value,s$name);want<-c(case_id='B8-U10-CENSORED-EVENTS',n='1600',target='B8-T03',horizon_ticks='12000',ticks_per_month='1000',groups='control|offered|nonrespondents',grid_ticks='0|3000|6000|9000|12000',ties='events_before_tied_censors')
need(setequal(names(spec),names(want))&&all(spec[names(want)]==want),'Unit10 R contract')
f<-read('data/anchor_participants.csv');b<-read('data/anchor_baseline.csv');old<-read('data/anchor_registry.csv');r<-read('data/anchor_observed_registry.csv')
for(d in list(f,b,old,r)) need(!anyDuplicated(d$student_id),'duplicate identifier')
need(nrow(f)==1600&&setequal(r$student_id,f$student_id),'invitation registry join');f<-f[order(f$student_id),];r<-r[match(f$student_id,r$student_id),]
need(all(f$responded_baseline %in% c('0','1')),'response flag');ids<-f$student_id[f$responded_baseline=='1'];need(setequal(b$student_id,ids)&&setequal(old$student_id,ids),'respondent join')
r$ticks<-vapply(r$exit_month,ticks,integer(1));r$group<-''
for(i in seq_len(nrow(r))) {
 x<-r[i,];need(x$cohort=='anchor'&&x$responded_baseline==f$responded_baseline[i]&&x$offer %in% c('0','1'),'cohort/response/offer')
 e<-x$exit_type;need(e %in% c('none','withdrawal','transfer','consent_withdrawn'),'event type');need(e!='none'||x$ticks==12000,'early no event');need(e!='consent_withdrawn'||x$ticks<12000,'consent horizon')
 known<-if(e=='consent_withdrawn') '0' else '1';label<-if(known=='0') '' else if(e=='withdrawal') '1' else '0';need(x$outcome_known==known&&x$withdrawal_label==label,'unknown/event label')
 if(x$responded_baseline=='1') {
  j<-match(x$student_id,b$student_id);k<-match(x$student_id,old$student_id)
  need(x$offer==b$offer[j]&&b$college_id[j]==f$college_id[i]&&x$exit_type==old$exit_type[k]&&x$exit_month==old$exit_month[k],'original history/assignment changed')
  need(x$observation_source=='v1.0.0_respondent_registry','respondent source');r$group[i]<-if(x$offer=='1') 'offered' else 'control'
 } else {
  need(x$offer=='0'&&e!='consent_withdrawn'&&x$observation_source=='v1.1.0_admin_no_offer','administrative policy');r$group[i]<-'nonrespondents'
 }
}
curve<-function(d,horizon=12000) {
 need(nrow(d)>0,'empty group');need(all(d$exit_type!='none'|d$ticks==12000)&&all(d$exit_type!='consent_withdrawn'|d$ticks<12000),'invalid censor horizon');times<-sort(unique(d$ticks[d$ticks<=horizon]));prob<-matrix(c(1,0,0),nrow=1);km<-1;rows<-list()
 for(t in times) {
  n<-sum(d$ticks>=t);z<-d[d$ticks==t,];dw<-sum(z$exit_type=='withdrawal');dt<-sum(z$exit_type=='transfer');dc<-sum(z$exit_type=='consent_withdrawn');da<-sum(z$exit_type=='none');before<-prob[1]
  trans<-diag(3);trans[1,]<-c(1-(dw+dt)/n,dw/n,dt/n);prob<-prob%*%trans;km<-km*(1-dw/n)
  need(abs(sum(prob)-1)<1e-12,'probability mass')
  rows[[length(rows)+1]]<-data.frame(ticks=t,risk=n,withdrawals=dw,transfers=dt,consent=dc,administrative=da,s_before=before,increment_w=before*dw/n,increment_t=before*dt/n,withdrawal=prob[2],transfer=prob[3],event_free=prob[1],withdrawal_only=1-km)
 }
 need(max(d$ticks)>=horizon||prob[1]==0,'follow-up does not support horizon')
 do.call(rbind,rows)
}
at<-function(a,t) {z<-c(withdrawal=0,transfer=0,event_free=1,withdrawal_only=0);ix<-which(a$ticks<=t);if(length(ix)) z[]<-unlist(a[max(ix),names(z)]);z}
o<-list();put<-function(k,v) o[[k]]<<-as.numeric(v);curves<-list();riskrows<-list();groups<-c('control','offered','nonrespondents');grid<-c(0,3000,6000,9000,12000)
for(g in c('pooled',groups)) {
 d<-if(g=='pooled') r else r[r$group==g,];n<-nrow(d);a<-curve(d);curves[[g]]<-a;w<-sum(d$exit_type=='withdrawal');tr<-sum(d$exit_type=='transfer');u<-sum(d$exit_type=='consent_withdrawn');admin<-n-w-tr-u
 values<-c(n=n,share=n/1600,withdrawals=w,transfers=tr,unknown=u,administrative=admin,lower=w/n,upper=(w+u)/n,risk_rows=nrow(a))
 for(k in names(values)) put(paste(g,k,sep='/'),values[k])
 for(t in grid) {v<-at(a,t);for(k in names(v)) put(paste(g,'month',t/1000,k,sep='/'),v[k])}
 riskrows[[g]]<-data.frame(group=g,a)
}
for(t in grid) for(k in c('withdrawal','transfer','event_free')) put(paste('standardized/month',t/1000,k,sep='/'),sum(vapply(groups,function(g) sum(r$group==g)/1600*at(curves[[g]],t)[k],numeric(1))))
put('diagnostic/pooled_minus_standardized',o[['pooled/month/12/withdrawal']]-o[['standardized/month/12/withdrawal']]);put('diagnostic/withdrawal_only_minus_competing',o[['pooled/month/12/withdrawal_only']]-o[['pooled/month/12/withdrawal']])
toy<-read('model-specs/unit10_histories.csv');need(nrow(toy)==14&&!anyDuplicated(paste(toy$fixture,toy$id))&&all(toy$fixture %in% c('eight','tied')),'toy identity');toy$ticks<-vapply(toy$exit_month,ticks,integer(1));need(all(toy$exit_type %in% c('none','withdrawal','transfer','consent_withdrawn')),'toy event');need(sum(toy$fixture=='eight')==8&&sum(toy$fixture=='tied')==6,'toy counts')
for(name in c('eight','tied')) {
 d<-toy[toy$fixture==name,];a<-curve(d);v<-at(a,12000);for(k in names(v)) put(paste('toy',name,k,sep='/'),v[k]);put(paste('toy',name,'lower',sep='/'),sum(d$exit_type=='withdrawal')/nrow(d));put(paste('toy',name,'upper',sep='/'),sum(d$exit_type %in% c('withdrawal','consent_withdrawn'))/nrow(d))
 if(name=='tied') for(j in seq_len(nrow(a))) for(k in names(a)) put(paste('toy/tied/step',j,k,sep='/'),a[j,k])
}
eight<-toy[toy$fixture=='eight',];est<-numeric(8)
for(j in seq_len(8)) {d<-eight;d$ticks[j]<-1000;d$exit_type[j]<-'consent_withdrawn';est[j]<-at(curve(d),12000)['withdrawal'];put(paste0('toy/uniform_mask/',d$id[j]),est[j])}
put('toy/uniform_mean',mean(est));put('toy/uniform_variance',mean((est-mean(est))^2));d<-eight;ix<-d$exit_type=='withdrawal';d$ticks[ix]<-1000;d$exit_type[ix]<-'consent_withdrawn';put('toy/selective_withdrawal',at(curve(d),12000)['withdrawal'])
keys<-sort(names(o));write.csv(data.frame(key=keys,value=vapply(o[keys],function(v)sprintf('%.17g',v),character(1))),output,row.names=FALSE,quote=FALSE,na='')
people<-r[,c('student_id','group','exit_month','ticks','exit_type','withdrawal_label','outcome_known','responded_baseline','offer','observation_source')]
write.csv(people,sub('\\.csv$','_people.csv',output),row.names=FALSE,quote=TRUE,na='')
risk<-do.call(rbind,riskrows);for(k in names(risk)[-1]) risk[[k]]<-vapply(risk[[k]],function(v)sprintf('%.17g',v),character(1));write.csv(risk,sub('\\.csv$','_risksets.csv',output),row.names=FALSE,quote=FALSE,na='');cat('UNIT10_BASE_R_OK\n')
