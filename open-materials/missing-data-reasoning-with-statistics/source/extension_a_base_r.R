# Independent beta-function predictive probabilities and native R summaries.
options(digits=17)
a<-commandArgs(trailingOnly=TRUE);stopifnot(length(a)==2);root<-a[1];out<-a[2]
readc<-function(n)read.csv(file.path(root,n),colClasses="character",na.strings=character(),check.names=FALSE)
cfg<-readc("model-specs/extension_a_r.csv");stopifnot(!anyDuplicated(cfg$name));s<-setNames(cfg$value,cfg$name)
stopifnot(s['n']=='6',s['low']=='1',s['high']=='5',s['prior_high']=='1',s['prior_low']=='1',s['replicates']=='32',s['max_m']=='800',s['grid']=='5|20|50|200|565|800',as.numeric(s['epsilon'])==.02,as.numeric(s['boundary_margin'])==1e-12)
r<-readc('data/extension_a.csv');stopifnot(identical(r$person,paste0('A',1:6)),all(r$score %in% c('','1','5')))
y<-as.numeric(r$score[r$score!='']);d<-6-length(y);v<-c(n=6,observed_n=length(y),missing_n=d,epsilon=.02);batches<-data.frame(replicate=numeric(),m=numeric(),mean=numeric(),B=numeric(),estimated_mcse=numeric(),exact_mcse=numeric(),error=numeric());draws<-data.frame(replicate=numeric(),draw=numeric(),K=numeric(),Q=numeric())
if(d==0){v<-c(v,posterior_mean=mean(y),posterior_variance=0,posterior_sd=0,required_m=0)}else{
 alpha<-1+sum(y==5);beta<-1+sum(y==1);k<-0:d;p<-choose(d,k)*exp(lbeta(alpha+k,beta+d-k)-lbeta(alpha,beta));q<-(sum(y)+d+4*k)/6;mu<-sum(p*q);vv<-sum(p*(q-mu)^2)
 v<-c(v,posterior_mean=mu,posterior_variance=vv,posterior_sd=sqrt(vv),required_m=max(2,ceiling(vv/.02^2)))
 for(i in seq_along(k)){v[paste0('law/',k[i],'/prob')]<-p[i];v[paste0('law/',k[i],'/q')]<-q[i]}
 t<-readc('random-inputs/extension_a_uniforms.csv');stopifnot(identical(as.integer(t$replicate),rep(1:32,each=800)),identical(as.integer(t$draw),rep(1:800,32)))
 u<-as.numeric(t$u);cdf<-cumsum(p);stopifnot(all(is.finite(u)),all(u>0 & u<1),all(vapply(cdf[-length(cdf)],function(c)all(abs(u-c)>1e-12),logical(1))))
 kk<-vapply(u,function(x)which(x<cdf)[1]-1L,integer(1));Q<-matrix(q[kk+1],nrow=32,byrow=TRUE);draws<-data.frame(replicate=as.integer(t$replicate),draw=as.integer(t$draw),K=kk,Q=q[kk+1])
 for(m in as.integer(strsplit(s['grid'],'|',fixed=TRUE)[[1]])){
  means<-apply(Q[,1:m,drop=FALSE],1,mean);between<-apply(Q[,1:m,drop=FALSE],1,var);errors<-means-mu;mc<-sqrt(between/m)
  batches<-rbind(batches,data.frame(replicate=1:32,m=m,mean=means,B=between,estimated_mcse=mc,exact_mcse=sqrt(vv/m),error=errors))
  v[paste0('grid/',m,'/replicate1_mean')]<-means[1];v[paste0('grid/',m,'/replicate1_error')]<-errors[1];v[paste0('grid/',m,'/estimated_mcse')]<-mc[1];v[paste0('grid/',m,'/exact_mcse')]<-sqrt(vv/m);v[paste0('grid/',m,'/empirical_rmse')]<-sqrt(mean(errors^2))
 }
}
write.csv(data.frame(key=sort(names(v)),value=unname(v[sort(names(v))])),out,row.names=FALSE,quote=FALSE)
write.csv(batches,sub('\\.csv$','_batches.csv',out),row.names=FALSE,quote=FALSE);write.csv(draws,sub('\\.csv$','_draws.csv',out),row.names=FALSE,quote=FALSE)
cat('EXTENSION_A_BASE_R_OK\n')
