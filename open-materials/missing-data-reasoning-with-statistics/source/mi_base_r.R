# Hand-checkable Rubin fixture is independent of the imputation case.
fixture_Q<-c(1,2,3);fixture_U<-c(4,4,4)
stopifnot(mean(fixture_U)==4,var(fixture_Q)==1,abs(mean(fixture_U)+(1+1/3)*var(fixture_Q)-16/3)<1e-14)
args <- commandArgs(trailingOnly=TRUE); root <- args[1]; out <- args[2]
spec <- read.csv(file.path(root,"model-specs/mi_fixture_r.csv"), colClasses="character")
stopifnot(identical(spec$name,c("model","observed_n","total_n","imputations","variance_df","prior")),identical(spec$value,c("iid_normal_intercept_only","6","8","8","5","1/sigma2")))
d <- read.csv(file.path(root,"data/mi_fixture.csv"), colClasses=c("character","numeric"), na.strings="")
stopifnot(identical(d$row_id,sprintf("F%02d",1:8)),sum(is.na(d$y))==2,identical(which(is.na(d$y)),7:8),all(is.finite(d$y[!is.na(d$y)])))
z <- read.csv(file.path(root,"random-inputs/mi_fixture_normals.csv"),colClasses="character")
stopifnot(nrow(z)==64,identical(z$imputation,as.character(rep(1:8,each=8))),all(z$chain=="1"),all(z$iteration=="1"),all(z$variable=="y"),all(z$variate_type=="standard_normal"),identical(z$update,rep(c(rep("variance",5),"mean","residual","residual"),8)),identical(z$row_key,rep(c(rep("parameter",6),"F07","F08"),8)),identical(z$position,rep(c(as.character(1:5),"1","1","1"),8)))
z$value <- as.numeric(z$value);stopifnot(all(is.finite(z$value)))
a <- mean(d$y,na.rm=TRUE); rss <- sum((d$y[!is.na(d$y)]-a)^2);stopifnot(rss>0)
r <- c(observed_mean=a,observed_rss=rss); Q <- U <- numeric(8)
for (m in 1:8) {
  s <- z$value[((m-1)*8+1):(m*8)];chi <- sum(s[1:5]^2);stopifnot(chi>0)
  s2 <- rss/chi;mu <- a+sqrt(s2/6)*s[6]; y <- d$y
  y[is.na(y)] <- mu+sqrt(s2)*s[7:8];Q[m] <- mean(y);U[m] <- var(y)/8
  p <- paste0("imputation/",m,"/")
  add <- c(chi2=chi,sigma2=s2,mu=mu);names(add)<-paste0(p,names(add));r<-c(r,add)
  add <- y;names(add)<-paste0(p,"completed/",d$row_id);r<-c(r,add)
  add <- c(Q=Q[m],U=U[m]);names(add)<-paste0(p,names(add));r<-c(r,add)
}
W <- mean(U);B <- var(Q);e <- (1+1/8)*B;T <- W+e
pooled<-c(Qbar=mean(Q),W=W,B=B,T=T,r=e/W,lambda=e/T,old_df=if(e==0) -1 else 7*(1+W/e)^2,MCSE_Qbar=sqrt(B/8));names(pooled)<-paste0("pool/",names(pooled));r<-c(r,pooled)
stopifnot(all(is.finite(r)),!anyDuplicated(names(r)))
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
write.table(data.frame(key=names(r),value=sprintf("%.17g",r)),out,sep=",",row.names=FALSE,col.names=TRUE,quote=FALSE)
cat("MI_BASE_R_OK end_to_end=true\n")
