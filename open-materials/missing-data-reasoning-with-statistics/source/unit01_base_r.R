# Independent base-R reconstruction from the invitation frame and baseline items.
a<-commandArgs(trailingOnly=TRUE);root<-a[1];outpath<-a[2]
read<-function(n)read.csv(file.path(root,'data',n),colClasses='character',na.strings=NULL,check.names=FALSE)
p<-read('anchor_participants.csv');p<-p[order(p$student_id),];b<-read('anchor_baseline.csv');items<-paste0('st',1:5)
stopifnot(nrow(p)==1600,length(unique(p$student_id))==1600,length(unique(b$student_id))==nrow(b),setequal(b$student_id,p$student_id[p$responded_baseline=='1']),all(p$responded_baseline %in% c('0','1')),all(p$part_time %in% c('0','1')))
idx<-match(p$student_id,b$student_id);stopifnot(all(p$college_id[!is.na(idx)]==b$college_id[idx[!is.na(idx)]]))
raw<-matrix('',nrow(p),5);raw[!is.na(idx),]<-as.matrix(b[idx[!is.na(idx)],items]);stopifnot(all(raw %in% c('','1','2','3','4','5')))
counts<-rowSums(raw!='');sums<-numeric(nrow(p));patterns<-character(nrow(p))
for(i in seq_len(nrow(p))){
 observed<-raw[i,]!='';values<-numeric(5)
 for(j in 1:5)if(observed[j])values[j]<-if(j==4)6-as.numeric(raw[i,j]) else as.numeric(raw[i,j])
 sums[i]<-sum(values);patterns[i]<-paste(as.integer(observed),collapse='')
}
n<-nrow(p);nr<-nrow(b);full<-counts==5;nc<-sum(full);known<-sum(counts);total<-sum(sums);cs<-sum(sums[full]);missing<-5*n-known
v<-c(invited=n,responded=nr,nonresponded=n-nr,complete_people=nc,partial_respondents=nr-nc,missing_full_scores=n-nc,response_fraction=nr/n,complete_fraction=nc/n,known_cells=known,missing_cells=missing,respondent_missing_cells=5*nr-known,nonrespondent_missing_cells=5*(n-nr),keyed_sum=total,complete_keyed_sum=cs,complete_mean=cs/(5*nc),item_lower=(total+missing)/(5*n),item_upper=(total+5*missing)/(5*n),person_lower=(cs+5*(n-nc))/(5*n),person_upper=(cs+25*(n-nc))/(5*n))
for(mask in 0:31){bits<-paste(as.integer(intToBits(mask))[5:1],collapse='');v[paste0('pattern/',bits)]<-sum(patterns==bits)}
for(g in c('0','1')){
 keep<-p$part_time==g;gn<-sum(keep);gr<-sum(p$responded_baseline[keep]=='1');x<-c(n=gn,responded=gr,complete=sum(full&keep),response_fraction=gr/gn);names(x)<-paste0('part_time_',g,'/',names(x));v<-c(v,x)
}
i<-which(full & raw[,4]!='3')[1];stopifnot(p$student_id[i]=='S1000118')
for(j in 1:5)v[paste0('example/',items[j])]<-as.numeric(raw[i,j])
v['example/keyed_st4']<-6-as.numeric(raw[i,4]);v['example/keyed_sum']<-sums[i];v['example/score']<-sums[i]/5;v['example/unkeyed_mean']<-sum(as.numeric(raw[i,]))/5

ex<-read.csv(file.path(root,'model-specs/unit01_exercises.csv'),colClasses='character',na.strings=NULL)
z<-as.character(ex[ex$record_id=='incomplete',paste0('x',1:5)]);present<-z!='';x<-numeric(5)
for(j in 1:5)if(present[j])x[j]<-if(j==4)6-as.numeric(z[j]) else as.numeric(z[j])
ns<-sum(present);ss<-sum(x);v['exercise/partial_keyed_sum']<-ss;v['exercise/available_mean']<-ss/ns;v['exercise/partial_lower']<-(ss+5-ns)/5;v['exercise/partial_upper']<-(ss+5*(5-ns))/5
ft<-ex[ex$record_id=='full_time',];pt<-ex[ex$record_id=='part_time',];v['exercise/full_time_rate']<-as.numeric(ft$x2)/as.numeric(ft$x1);v['exercise/part_time_rate']<-as.numeric(pt$x2)/as.numeric(pt$x1);v['exercise/respondent_part_time_share']<-as.numeric(pt$x2)/(as.numeric(pt$x2)+as.numeric(ft$x2))
v<-v[order(names(v))];writeLines(c('key,value',paste(names(v),sprintf('%.17g',v),sep=',')),outpath)
people<-c('student_id,college_id,responded_baseline,observed_items,keyed_sum,complete_score,pattern',vapply(seq_len(n),function(i)paste(p$student_id[i],p$college_id[i],p$responded_baseline[i],counts[i],sums[i],if(full[i])sprintf('%.17g',sums[i]/5) else '',patterns[i],sep=','),character(1)))
writeLines(people,sub('\\.csv$','_people.csv',outpath));cat('UNIT01_BASE_R_OK\n')
