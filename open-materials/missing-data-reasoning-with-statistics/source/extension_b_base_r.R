# Independent sequential posterior-predictive recursion, preserving raw items.
options(digits=17)
a<-commandArgs(trailingOnly=TRUE);stopifnot(length(a)==2);root<-a[1];out<-a[2]
readc<-function(n)read.csv(file.path(root,n),colClasses='character',na.strings=character(),check.names=FALSE)
cfg<-readc('model-specs/extension_b_r.csv');stopifnot(!anyDuplicated(cfg$name));s<-setNames(cfg$value,cfg$name)
stopifnot(s['n']=='6',s['prior_per_profile']=='1',s['reverse_item']=='st4',s['profiles']=='A|B|C|D',s['items']=='st1|st2|st3|st4|st5')
items<-paste0('st',1:5);r<-readc('data/extension_b.csv');pr<-readc('model-specs/extension_b_profiles.csv');stopifnot(identical(r$person,paste0('B',1:6)),identical(pr$profile,LETTERS[1:4]),all(unlist(r[items]) %in% c('',as.character(1:5))),all(unlist(pr[items]) %in% as.character(1:5)))
raw<-matrix(as.numeric(as.matrix(pr[items])),nrow=4);stopifnot(nrow(unique(raw))==4);key<-raw;key[,4]<-6-key[,4];scores<-rowMeans(key);allowed<-vector('list',6);nobs<-integer(6)
for(i in 1:6){nobs[i]<-sum(r[i,items]!='');allowed[[i]]<-which(vapply(1:4,function(k)all(vapply(1:5,function(j)r[i,items[j]]=='' || as.numeric(r[i,items[j]])==raw[k,j],logical(1))),logical(1)));stopifnot(length(allowed[[i]])>0)}
known<-which(nobs==5);unknown<-which(nobs<5);alpha<-rep(1,4);z<-integer(6)
for(i in known){z[i]<-allowed[[i]][1];alpha[z[i]]<-alpha[z[i]]+1}
assignments<-list();masses<-numeric();ok<-logical()
walk<-function(j,counts,prob,z){
 if(j>length(unknown)){assignments[[length(assignments)+1]]<<-z;masses[length(masses)+1]<<-prob;ok[length(ok)+1]<<-all(vapply(unknown,function(i)z[i] %in% allowed[[i]],logical(1)));return()}
 i<-unknown[j]
 for(k in 1:4){new<-counts;new[k]<-new[k]+1;zz<-z;zz[i]<-k;walk(j+1,new,prob*counts[k]/sum(counts),zz)}
}
walk(1,alpha,1,z);zz<-do.call(rbind,assignments);Q<-apply(zz,1,function(x)mean(scores[x]));evidence<-sum(masses[ok]);pi<-ifelse(ok,masses/evidence,0);ps<-masses/sum(masses)
v<-c('data/n'=6,'data/full_scores'=length(known),'data/partial_people'=sum(nobs>0 & nobs<5),'data/observed_item_cells'=sum(nobs),'data/discarded_partial_cells'=sum(nobs[unknown]),'data/state_count'=length(Q))
for(route in c('item','score')){p<-if(route=='item')pi else ps;mu<-sum(p*Q);v[paste0(route,'/evidence')]<-if(route=='item')evidence else sum(masses);v[paste0(route,'/mean')]<-mu;v[paste0(route,'/variance')]<-sum(p*(Q-mu)^2)}
people<-data.frame();lo<-hi<-available<-numeric()
for(i in 1:6){vals<-vapply(which(r[i,items]!=''),function(j)if(j==4)6-as.numeric(r[i,items[j]])else as.numeric(r[i,items[j]]),numeric(1));lo[i]<-(sum(vals)+5-length(vals))/5;hi[i]<-(sum(vals)+5*(5-length(vals)))/5;if(length(vals))available<-c(available,mean(vals));row<-data.frame(person=r$person[i],observed_items=nobs[i],logical_lower=lo[i],logical_upper=hi[i])
 for(route in c('item','score')){p<-if(route=='item')pi else ps;mu<-sum(p*scores[zz[,i]]);row[[paste0(route,'_mean')]]<-mu;row[[paste0(route,'_variance')]]<-sum(p*(scores[zz[,i]]-mu)^2)};people<-rbind(people,row)}
v<-c(v,logical_lower=mean(lo),logical_upper=mean(hi),available_mean_n=length(available),available_mean=mean(available))
states<-data.frame(assignment=apply(zz[,unknown,drop=FALSE],1,function(x)paste0(LETTERS[x],collapse='')),Q=Q,item_probability=pi,score_probability=ps)
write.csv(data.frame(key=sort(names(v)),value=unname(v[sort(names(v))])),out,row.names=FALSE,quote=FALSE);write.csv(people,sub('\\.csv$','_people.csv',out),row.names=FALSE,quote=FALSE);write.csv(states,sub('\\.csv$','_states.csv',out),row.names=FALSE,quote=FALSE);cat('EXTENSION_B_BASE_R_OK\n')
