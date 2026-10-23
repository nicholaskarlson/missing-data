a<-commandArgs(trailingOnly=TRUE);r<-read.csv(file.path(a[1],'data/anchor_observed_registry.csv'),colClasses='character',na.strings=NULL)
stopifnot(nrow(r)==1600,length(unique(r$student_id))==1600,all(r$outcome_known %in% c('0','1')),all(r$withdrawal_label[r$outcome_known=='0']==''),all(r$withdrawal_label[r$outcome_known=='1'] %in% c('0','1')))
groups<-list(all=rep(TRUE,nrow(r)),control=r$responded_baseline=='1' & r$offer=='0',offered=r$responded_baseline=='1' & r$offer=='1',nonrespondents=r$responded_baseline=='0');out<-numeric()
for(k in names(groups)){
 g<-r[groups[[k]],];n<-nrow(g);known<-sum(g$outcome_known=='1');events<-sum(g$withdrawal_label=='1');unknown<-n-known
 values<-c(n=n,known=known,unknown=unknown,events=events,known_label_fraction=events/known,lower=events/n,upper=(events+unknown)/n)
 names(values)<-paste0(k,'/',names(values));out<-c(out,values)
}
out['diagnostic/response_share_standardization']<-(1219/1600)*out['control/known_label_fraction']+(381/1600)*out['nonrespondents/known_label_fraction']
out<-out[order(names(out))];writeLines(c('key,value',paste(names(out),sprintf('%.17g',out),sep=',')),a[2]);cat('WITHDRAWAL_BASE_R_OK\n')
