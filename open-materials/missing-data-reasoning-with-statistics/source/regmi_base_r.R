# Independent base-R normal-regression MI; lm.fit/chol/qt, no add-on packages.
a<-commandArgs(trailingOnly=TRUE);root<-a[1];output<-a[2]
cfg<-read.csv(file.path(root,'model-specs/regmi_r.csv'),colClasses='character');stopifnot(!anyDuplicated(cfg$name));cfg<-setNames(cfg$value,cfg$name);n<-as.integer(cfg[['n']]);m<-as.integer(cfg[['imputations']]);stopifnot(cfg[['p']]=='2',m>=2,cfg[['prior']]=='1/sigma2',cfg[['target']]=='conditional_model_slope',cfg[['support']]=='unbounded',cfg[['gram_relative_floor']]=='1e-12')
raw<-read.csv(file.path(root,'data/regmi_calibration.csv'),colClasses='character',na.strings=NULL);raw<-raw[order(raw$row_id),];stopifnot(identical(raw$row_id,sprintf('P%03d',seq_len(n))))
x<-as.numeric(raw$x);seen<-raw$y!='';y<-rep(NA_real_,n);y[seen]<-as.numeric(raw$y[seen]);stopifnot(all(is.finite(x)),all(is.finite(y[seen])))
fit<-function(x,y){
 D<-cbind(1,x);g<-crossprod(D);stopifnot(length(y)>2,all(is.finite(y)),all(is.finite(D)),det(g)>1e-12*g[1,1]*g[2,2])
 f<-lm.fit(D,y);stopifnot(f$rank==2);ss<-sum(f$residuals^2);stopifnot(is.finite(ss),ss>1e-24*max(1,sum(y*y)))
 list(beta=f$coefficients,rss=ss,V=solve(g))
}
pool<-function(q,u,dfcom){
 stopifnot(length(q)>=2,length(q)==length(u),all(is.finite(q)),all(is.finite(u)),all(u>0),dfcom>0)
 m<-length(q);meanq<-mean(q);w<-mean(u);b<-var(q);extra<-(1+1/m)*b;total<-w+extra;relative<-extra/w;lambda<-extra/total;old<-if(b>0)(m-1)/lambda^2 else -1;obs<-dfcom*(dfcom+1)/(dfcom+3)*(1-lambda);df<-1/((if(old>0)1/old else 0)+1/obs);stopifnot(is.finite(df),df>=1)
 se<-sqrt(total);t<-qt(.975,df);c(qbar=meanq,W=w,B=b,T=total,SE=se,r=relative,lambda=lambda,dfcom=dfcom,dfold=old,dfobs=obs,df=df,fmi=(relative+2/(df+3))/(relative+1),MCSE=sqrt(b/m),critical=t,lower=meanq-t*se,upper=meanq+t*se)
}
f<-fit(x[seen],y[seen]);d<-sum(seen)-2;v<-c('input/n'=n,'input/nobs'=sum(seen),'input/nmis'=sum(!seen),'fit/beta0'=f$beta[1],'fit/beta1'=f$beta[2],'fit/RSS'=f$rss,'fit/df'=d,'fit/V00'=f$V[1,1],'fit/V01'=f$V[1,2],'fit/V11'=f$V[2,2],'complete_case/Q'=f$beta[2],'complete_case/U'=f$rss/d*f$V[2,2])
# Avoid inherited lm coefficient names in keyed output.
names(v)<-c('input/n','input/nobs','input/nmis','fit/beta0','fit/beta1','fit/RSS','fit/df','fit/V00','fit/V01','fit/V11','complete_case/Q','complete_case/U')
if(all(seen)){
 q<-unname(f$beta[2]);u<-f$rss/d*f$V[2,2];t<-qt(.975,d);p<-c(qbar=q,W=u,B=0,T=u,SE=sqrt(u),r=0,lambda=0,dfcom=d,dfold=-1,dfobs=d,df=d,fmi=0,MCSE=0,critical=t,lower=q-t*sqrt(u),upper=q+t*sqrt(u));v<-c(v,setNames(p,paste0('pool/',names(p))),'branch/no_missing'=1)
}else{
 tape<-read.csv(file.path(root,'random-inputs/regmi_normals.csv'),colClasses='character');wanted<-character(0)
 for(j in seq_len(m))for(role in c('variance','beta','residual')){
  idx<-switch(role,variance=as.character(seq_len(d)),beta=c('0','1'),residual=raw$row_id[!seen]);wanted<-c(wanted,paste(j,role,idx,sep='|'))
 }
 stopifnot(identical(paste(tape$imputation,tape$role,tape$index,sep='|'),wanted),all(is.finite(as.numeric(tape$value))))
 tape$value<-as.numeric(tape$value);l<-t(chol(f$V));q<-u<-numeric(m)
 for(j in seq_len(m)){
  z<-tape[tape$imputation==as.character(j),];variancez<-z$value[z$role=='variance'];chi<-sum(variancez^2);stopifnot(chi>0);s2<-f$rss/chi;beta<-f$beta+sqrt(s2)*as.vector(l%*%z$value[z$role=='beta']);filled<-y;filled[!seen]<-beta[1]+beta[2]*x[!seen]+sqrt(s2)*z$value[z$role=='residual'];g<-fit(x,filled);q[j]<-g$beta[2];u[j]<-g$rss/(n-2)*g$V[2,2];prefix<-paste0('imputation/',j,'/');p<-c(chi,s2,beta[1],beta[2],q[j],u[j],g$rss);names(p)<-paste0(prefix,c('chi','sigma2','beta0','beta1','Q','U','RSS'));v<-c(v,p,setNames(filled,paste0(prefix,'completed/',raw$row_id)))
 }
 filled<-y;filled[!seen]<-f$beta[1]+f$beta[2]*x[!seen];g<-fit(x,filled);v['single/Q']<-unname(g$beta[2]);v['single/U']<-g$rss/(n-2)*g$V[2,2];p<-pool(q,u,n-2);v<-c(v,setNames(p,paste0('pool/',names(p))),'branch/no_missing'=0)
}
stopifnot(!anyDuplicated(names(v)),all(is.finite(v)));v<-v[order(names(v))];writeLines(c('key,value',paste(names(v),sprintf('%.17g',v),sep=',')),output);cat('REGMI_BASE_R_OK calibration_only=true\n')
