# Independent bounded item completion and finite-roster posterior moments, base R only.
a<-commandArgs(trailingOnly=TRUE);root<-a[1];path<-a[2]
spec<-read.csv(file.path(root,'model-specs/unit05_r.csv'),colClasses='character');stopifnot(!anyDuplicated(spec$name));spec<-setNames(spec$value,spec$name);n<-as.integer(spec[['n']]);ndraw<-as.integer(spec[['imputations']]);stopifnot(spec[['prior_per_category']]=='1',ndraw>=2,spec[['target']]=='finite_roster_mean',spec[['groups']]=='00|01|10|11',spec[['boundary_margin']]=='1e-10',spec[['dependence_flag']]=='.10',spec[['reverse_item']]=='st4')
read<-function(name)read.csv(file.path(root,'data',name),colClasses='character',na.strings=NULL)
p<-read('anchor_participants.csv');p<-p[order(p$student_id),];b<-read('anchor_baseline.csv');stopifnot(nrow(p)==n,!anyDuplicated(p$student_id),!anyDuplicated(b$student_id),setequal(b$student_id,p$student_id[p$responded_baseline=='1']))
stopifnot(all(p$part_time %in% c('0','1')),all(p$first_gen %in% c('0','1')),all(p$responded_baseline %in% c('0','1')))
items<-paste0('st',1:5);levels<-c('00','01','10','11');group<-paste0(p$part_time,p$first_gen);stopifnot(setequal(group,levels));idx<-match(p$student_id,b$student_id);stopifnot(all(p$college_id[!is.na(idx)]==b$college_id[idx[!is.na(idx)]]))
raw<-matrix('',n,5);raw[!is.na(idx),]<-as.matrix(b[idx[!is.na(idx)],items]);stopifnot(all(raw %in% c('','1','2','3','4','5')));seen<-raw!='';keyed<-matrix(0,n,5);keyed[seen]<-as.numeric(raw[seen]);keyed[,4]<-ifelse(seen[,4],6-keyed[,4],0);miss<-which(!seen,arr.ind=TRUE);miss<-miss[order(miss[,1],miss[,2]),,drop=FALSE];full<-rowSums(seen)==5;den<-5*n;known<-sum(seen);s<-sum(keyed)
count<-array(0,c(4,5,5));nm<-matrix(0,4,5);alpha<-array(0,c(4,5,5));v<-c(invited=n,known_cells=known,missing_cells=nrow(miss),keyed_sum=s,complete_people=sum(full),'support/lower'=(s+nrow(miss))/den,'support/upper'=(s+5*nrow(miss))/den);expect<-s;residvar<-paramvar<-0
for(g in 1:4)for(j in 1:5){
 keep<-group==levels[g];count[g,j,]<-tabulate(keyed[keep&seen[,j],j],nbins=5);stopifnot(sum(count[g,j,])>0);nm[g,j]<-sum(keep&!seen[,j]);alpha[g,j,]<-count[g,j,]+1;A<-sum(alpha[g,j,]);mu<-sum((1:5)*alpha[g,j,])/A;vv<-sum((1:5)^2*alpha[g,j,])/A-mu^2;rr<-nm[g,j]*A*vv/(A+1);pp<-nm[g,j]^2*vv/(A+1);expect<-expect+nm[g,j]*mu;residvar<-residvar+rr;paramvar<-paramvar+pp;prefix<-paste0('cell/',levels[g],'/',items[j],'/')
 for(k in 1:5){v[paste0(prefix,'count/',k)]<-count[g,j,k];v[paste0(prefix,'alpha/',k)]<-alpha[g,j,k]}
 z<-c(A=A,missing=nm[g,j],mean=mu,variance=vv,residual_variance_total=rr,parameter_variance_total=pp);v<-c(v,setNames(z,paste0(prefix,names(z))))
}
mean<-expect/den;variance<-(residvar+paramvar)/den^2;stopifnot(is.finite(variance),variance>=0);sd<-sqrt(variance);v<-c(v,'posterior/mean'=mean,'posterior/variance'=variance,'posterior/SD'=sd,'posterior/residual_variance'=residvar/den^2,'posterior/parameter_variance'=paramvar/den^2,'posterior/chebyshev_lower'=max(v[['support/lower']],mean-sqrt(variance/.05)),'posterior/chebyshev_upper'=min(v[['support/upper']],mean+sqrt(variance/.05)))
center<-keyed[full,,drop=FALSE];fg<-group[full]
for(g in levels){stopifnot(sum(fg==g)>=2);center[fg==g,]<-sweep(center[fg==g,,drop=FALSE],2,colMeans(center[fg==g,,drop=FALSE]))}
cross<-crossprod(center);correlation<-cross/sqrt(outer(diag(cross),diag(cross)));stopifnot(all(is.finite(correlation)));corrs<-numeric(0)
for(j in 1:4)for(k in (j+1):5){v[paste0('diagnostic/correlation/',j,k)]<-correlation[j,k];corrs<-c(corrs,abs(correlation[j,k]))}
v['diagnostic/max_abs_correlation']<-max(corrs);v['diagnostic/dependence_flag']<-as.integer(max(corrs)>.10);people<-list()
if(nrow(miss)==0){v<-c(v,'branch/no_missing'=1,'MC/Qbar'=mean,'MC/B'=0,'MC/estimated_SE'=0,'MC/exact_SE'=0)}else{
 tape<-read.csv(file.path(root,'random-inputs/unit05_uniforms.csv'),colClasses='character');postkeys<-unlist(lapply(levels,function(g)unlist(lapply(items,function(j)paste(g,j,1:5,sep='/')))));cellkeys<-paste(p$student_id[miss[,1]],items[miss[,2]],sep='/');expected<-unlist(lapply(seq_len(ndraw),function(d)c(paste(d,'posterior',postkeys,sep='|'),paste(d,'cell',cellkeys,sep='|'))));stopifnot(identical(paste(tape$imputation,tape$role,tape$key,sep='|'),expected));tape$value<-as.numeric(tape$value);stopifnot(all(is.finite(tape$value)),all(tape$value>0&tape$value<1));Q<-numeric(ndraw)
 for(d in seq_len(ndraw)){
  z<-tape[tape$imputation==as.character(d),];u<-z$value[z$role=='posterior'];theta<-gamma<-array(0,c(4,5,5));cursor<-1
  for(g in 1:4)for(j in 1:5){
   gg<-qgamma(u[cursor:(cursor+4)],shape=alpha[g,j,],scale=1);stopifnot(all(is.finite(gg)),all(gg>0));cursor<-cursor+5;gamma[g,j,]<-gg;theta[g,j,]<-gg/sum(gg)
   for(k in 1:5){v[paste0('draw/',d,'/gamma/',levels[g],'/',items[j],'/',k)]<-gg[k];v[paste0('draw/',d,'/theta/',levels[g],'/',items[j],'/',k)]<-theta[g,j,k]}
  }
  filled<-keyed;uu<-z$value[z$role=='cell']
  for(h in seq_len(nrow(miss))){i<-miss[h,1];j<-miss[h,2];g<-match(group[i],levels);cum<-cumsum(theta[g,j,]);stopifnot(all(abs(uu[h]-cum[1:4])>1e-10));category<-which(uu[h]<cum)[1];stopifnot(!is.na(category));filled[i,j]<-category}
  total<-rowSums(filled);Q[d]<-sum(total)/den;v[paste0('draw/',d,'/Q')]<-Q[d];rawfilled<-filled;rawfilled[,4]<-6-rawfilled[,4];people[[d]]<-data.frame(imputation=as.character(d),student_id=p$student_id,group=group,st1=as.character(rawfilled[,1]),st2=as.character(rawfilled[,2]),st3=as.character(rawfilled[,3]),st4=as.character(rawfilled[,4]),st5=as.character(rawfilled[,5]),score=sprintf('%.17g',total/5),stringsAsFactors=FALSE)
  if(d==1){i<-match(spec[['worked_example_id']],p$student_id);stopifnot(!is.na(i));h<-which(miss[,1]==i)[1];v['example/uniform']<-uu[h];v['example/score']<-total[i]/5;for(j in 1:5)v[paste0('example/filled_',items[j])]<-rawfilled[i,j]}
 }
 v<-c(v,'branch/no_missing'=0,'MC/Qbar'=mean(Q),'MC/B'=var(Q),'MC/estimated_SE'=sqrt(var(Q)/ndraw),'MC/exact_SE'=sd/sqrt(ndraw))
}
stopifnot(!anyDuplicated(names(v)),all(is.finite(v)));v<-v[order(names(v))];writeLines(c('key,value',paste(names(v),sprintf('%.17g',v),sep=',')),path)
if(length(people))write.table(do.call(rbind,people),sub('\\.csv$','_people.csv',path),sep=',',row.names=FALSE,col.names=TRUE,quote=FALSE) else writeLines('imputation,student_id,group,st1,st2,st3,st4,st5,score',sub('\\.csv$','_people.csv',path))
cat('UNIT05_BASE_R_OK finite_roster=true\n')
