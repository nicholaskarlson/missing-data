# Independent base-R census rows and exhaustive form allocations, no analysis packages.
a<-commandArgs(trailingOnly=TRUE);root<-a[1];output<-a[2]
rd<-function(p)read.csv(file.path(root,p),colClasses='character',check.names=FALSE,na.strings=NULL)
need<-function(x,m)if(!isTRUE(x))stop(m)
items<-paste0('wb',1:5);audit<-c(items,paste0('st',1:5));c<-rd('model-specs/unit12_r.csv');need(!anyDuplicated(c$name),'duplicate specification');s<-setNames(c$value,c$name)
need(identical(as.character(s[c('n_randomized','n_respondents','wave','items','audit_items','reverse_items','forms','support','design','alpha','toy_n','fixtures')]),c('1219','956','6',paste(items,collapse='|'),paste(audit,collapse='|'),'wb4|st4','1:wb3,st3|2:wb4,st4|3:wb5,st5','1|5','independent_uniform_three_forms','.05','4','parallel|mixed|unequal_blocks')),'fixed Unit12 R specification')
f<-rd('data/anchor_participants.csv');b<-rd('data/anchor_baseline.csv');v<-rd('data/anchor_followup.csv')
need(nrow(f)==1600&&!anyDuplicated(f$student_id)&&all(f$responded_baseline%in%c('0','1')),'frame identity');need(nrow(b)==1219&&!anyDuplicated(b$student_id)&&setequal(b$student_id,f$student_id[f$responded_baseline=='1']),'randomized membership');need(all(b$offer%in%c('0','1'))&&all(b$college_id==f$college_id[match(b$student_id,f$student_id)]),'offer/college join');need(length(unique(b$college_id))==16,'college count')
for(g in unique(b$college_id)){nn<-table(factor(b$offer[b$college_id==g],levels=c('0','1')));need(all(nn>0)&&abs(nn[1]-nn[2])<=1,'original balance')}
keys<-paste(v$student_id,v$wave_month);need(!anyDuplicated(keys)&&setequal(keys,as.vector(outer(b$student_id,c('3','6','12'),paste))),'wave inventory');need(all(v$responded%in%c('0','1'))&&all(as.matrix(v[audit])%in%c('',as.character(1:5))),'category support')
for(i in seq_len(nrow(v))){r<-v[i,];if(r$responded=='0'){need(r$form=='0'&&all(r[audit]==''),'wave nonresponse answers/form')}else if(r$wave_month=='6'){need(r$form%in%c('1','2','3'),'unknown planned form');j<-as.integer(r$form)+2;need(r[[audit[j]]]==''&&r[[audit[j+5]]]=='','planned omission filled')}else need(r$form=='0','unexpected form')}
b<-b[order(b$student_id),];six<-v[v$wave_month=='6',];six<-six[match(b$student_id,six$student_id),];need(sum(six$responded=='1')==956,'fixed respondent roster');p<-data.frame(student_id=b$student_id,college_id=b$college_id,offer=b$offer,responded6=six$responded,form=six$form,six[audit],stringsAsFactors=FALSE)
T<-rd('model-specs/unit12_teaching.csv');need(nrow(T)==12&&!anyDuplicated(paste(T$fixture,T$id))&&setequal(paste(T$fixture,T$id),as.vector(outer(c('parallel','mixed','unequal_blocks'),paste0('P',1:4),paste)))&&all(as.matrix(T[items])%in%as.character(1:5)),'teaching inventory')
X<-matrix(as.numeric(as.matrix(p[audit])),ncol=10);X[,c(4,9)]<-6-X[,c(4,9)];known<-!is.na(X);respond<-p$responded6=='1';planned<-matrix(FALSE,nrow(p),10)
for(form in 1:3){planned[p$form==as.character(form),form+2]<-TRUE;planned[p$form==as.character(form),form+7]<-TRUE}
wave<-matrix(rep(!respond,10),ncol=10);incidental<-!known&!planned&!wave;weights<-matrix(rep(c(1,1,1.5,1.5,1.5),each=nrow(p)),ncol=5);weights[planned[,1:5]]<-0;total<-rowSums(X[,1:5]*weights,na.rm=TRUE);unknown<-rowSums(is.na(X[,1:5])*weights);Hlo<-(total+unknown)/5;Hhi<-(total+5*unknown)/5;ks<-rowSums(X[,1:5],na.rm=TRUE);mi<-5-rowSums(known[,1:5]);fullLo<-(ks+mi)/5;fullHi<-(ks+5*mi)/5
for(label in c('wb','st')){ix<-if(label=='wb')1:5 else 6:10;for(name in c('observed','planned','wave','incidental')){A<-switch(name,observed=known,planned=planned,wave=wave,incidental=incidental);p[[paste0(label,'_',name)]]<-rowSums(A[,ix])}}
p$keyed_wb_sum<-ks;p$full_score<-'';p$full_lower<-fullLo;p$full_upper<-fullHi;p$H_known_total<-ifelse(respond,total,NA);p$H_lower<-ifelse(respond,Hlo,NA);p$H_upper<-ifelse(respond,Hhi,NA)
cells<-list();k<-1
for(form in 0:3)for(j in 1:10){ix<-p$form==as.character(form);cells[[k]]<-data.frame(form=as.character(form),item=audit[j],roster_n=sum(ix),observed_n=sum(known[ix,j]),planned_n=sum(planned[ix,j]),incidental_n=sum(incidental[ix,j]),wave_n=sum(wave[ix,j]),keyed_sum=sum(X[ix,j],na.rm=TRUE));k<-k+1};cells<-do.call(rbind,cells)
o<-c(n_randomized=nrow(p),n_respondents=sum(respond),n_wave_missing=sum(!respond),cells=nrow(cells))
for(name in c('observed','planned','incidental','wave'))o<-c(o,setNames(sum(cells[[paste0(name,'_n')]]),paste0('audit/',name,'_n')))
for(form in 1:3)o<-c(o,setNames(sum(p$form==as.character(form)),paste0('form/',form,'/n')))
for(pop in c('respondents','randomized')){ix<-if(pop=='respondents')respond else rep(TRUE,nrow(p));o<-c(o,setNames(c(mean(fullLo[ix]),mean(fullHi[ix])),paste0('support/',pop,c('/lower','/upper'))))}
n<-sum(respond);vmax<-8/(25*n);radius<-sqrt(vmax/.05);o<-c(o,'H/lower'=mean(Hlo[respond]),'H/upper'=mean(Hhi[respond]),'design/variance_upper'=vmax,'design/radius'=radius,'design/lower'=max(1,mean(Hlo[respond])-radius),'design/upper'=min(5,mean(Hhi[respond])+radius),'burden/omitted_fraction'=.2)
F<-expand.grid(rep(list(1:3),4));code<-apply(F,1,paste0,collapse='');F<-as.matrix(F[order(code),]);masks<-list();index<-1
for(fixture in c('parallel','mixed','unequal_blocks')){tt<-T[T$fixture==fixture,];tt<-tt[order(tt$id),];Y<-matrix(as.numeric(as.matrix(tt[items])),ncol=5);Y[,4]<-6-Y[,4];target<-mean(Y);H<-numeric(81);V<-numeric(81);naive<-numeric(81)
 for(a in 1:81){h<-numeric(4);nv<-numeric(4);pairDiff<-numeric(4)
  for(i in 1:4){pair<-Y[i,(3:5)[(1:3)!=F[a,i]]];h[i]<-(sum(Y[i,1:2])+1.5*sum(pair))/5;nv[i]<-(sum(Y[i,1:2])+sum(pair))/4;pairDiff[i]<-(pair[1]-pair[2])^2}
  H[a]<-mean(h);naive[a]<-mean(nv);V[a]<-3*sum(pairDiff)/1600;masks[[index]]<-data.frame(fixture=fixture,forms=paste0(F[a,],collapse=''),target=target,H=H[a],naive=naive[a],Vhat=V[a]);index<-index+1
 }
 selective<-mean((rowSums(Y[,1:2])+1.5*(rowSums(Y[,3:5])-apply(Y[,3:5],1,max)))/5)
 values<-c(target=target,mean_H=mean(H),variance=mean((H-mean(H))^2),mean_Vhat=mean(V),mean_naive=mean(naive),selective_H=selective,assignments=81);o<-c(o,setNames(values,paste0('toy/',fixture,'/',names(values))))
}
triples<-expand.grid(rep(list(1:5),3));ps<-(triples[,1]-triples[,2])^2+(triples[,1]-triples[,3])^2+(triples[,2]-triples[,3])^2;o<-c(o,'bound/max_pair_square'=max(ps),'bound/triples'=nrow(triples));need(all(is.finite(o))&&!anyDuplicated(names(o)),'invalid output');options(digits=17)
write.table(data.frame(key=names(o)[order(names(o))],value=sprintf('%.17g',o[order(names(o))])),output,sep=',',row.names=FALSE,quote=FALSE)
for(name in c('people','cells','masks')){x<-switch(name,people=p,cells=cells,masks=do.call(rbind,masks));write.table(x,sub('\\.csv$',paste0('_',name,'.csv'),output),sep=',',row.names=FALSE,quote=FALSE,na='')}
cat('UNIT12_BASE_R_OK conditional_design=true\n')
