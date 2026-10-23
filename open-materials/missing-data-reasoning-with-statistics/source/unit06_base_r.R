# Independent base R: log-beta expansion, exact moments, joint draws and Gibbs matrices.
args <- commandArgs(trailingOnly=TRUE);root <- args[1];output <- args[2]
read <- function(n) read.csv(file.path(root,n),colClasses='character',check.names=FALSE,na.strings=NULL)
rs <- read('model-specs/unit06_r.csv');expected <- c(case_id='B8-U06-JOINT-BINARY-001',n='1600',prior_per_cell='1',imputations='20',high_threshold='4',groups='00|01|10|11',state_order='00|01|10|11',boundary_margin='1e-10',mixing_warning='1e-8',sweep_grid='1|2|5|10|25|50')
if(anyDuplicated(rs$name) || !identical(setNames(rs$value,rs$name),expected)) stop('fixed R contract')
states <- c('00','01','10','11');labels <- c('00','01','10','11');rolelabels <- c('row0','row1','col0','col1','both');grid <- c(1,2,5,10,25,50);margin <- 1e-10
out <- numeric();put <- function(k,v) {if(length(v)!=1 || !is.finite(v)) stop('nonfinite output');out[k] <<- v}
fmt <- function(v) sprintf('%.17g',v)
kernel <- function(a,b,order) {
 if(!order %in% c('AB','BA')) stop('scan order')
 if(length(a)!=2 || length(b)!=2 || any(!is.finite(c(a,b))) || any(c(a,b)<=0 | c(a,b)>=1)) stop('conditional support')
 t <- matrix(0,4,4)
 for(i in 0:3) for(j in 0:3) {pa <- a[(if(order=='AB') i%%2 else j%%2)+1];pb <- b[(if(order=='AB') j%/%2 else i%/%2)+1];t[i+1,j+1] <- (if(j%/%2==1) pa else 1-pa)*(if(j%%2==1) pb else 1-pb)}
 t
}
stationary <- function(t) {a <- t(t-diag(4));a[4,] <- 1;as.vector(solve(a,c(0,0,0,1)))}
audit <- function(prefix,p) {
 a <- c(p[3]/(p[1]+p[3]),p[4]/(p[2]+p[4]));b <- c(p[2]/(p[1]+p[2]),p[4]/(p[3]+p[4]))
 for(j in 1:4) put(paste0(prefix,'/',c('a0','a1','b0','b1')[j]),c(a,b)[j])
 for(order in c('AB','BA')) {
  t <- kernel(a,b,order);power <- diag(4)
  for(i in 1:4) for(j in 1:4) put(paste0(prefix,'/',order,'/transition/',i-1,'/',j-1),t[i,j])
  put(paste0(prefix,'/',order,'/invariant_residual'),max(abs(as.vector(p%*%t)-p)))
  for(sweep in 1:50) {power <- power%*%t;if(sweep %in% grid) put(paste0(prefix,'/',order,'/maxTV/',sweep),max(rowSums(abs(sweep(power,2,p,'-')))/2))}
 }
}
mixture <- function(co,pa,both) {
 if(length(co)!=4 || length(pa)!=4 || any(!is.finite(c(co,pa,both))) || any(c(co,pa,both)<0) || any(c(co,pa,both)!=as.integer(c(co,pa,both)))) stop('counts')
 z <- as.matrix(expand.grid(lapply(rev(pa),function(k) 0:k),KEEP.OUT.ATTRS=FALSE))[,4:1,drop=FALSE]
 a <- z[,1];b <- z[,2];c <- z[,3];d <- z[,4];alpha <- cbind(co[1]+1+a+c,co[2]+1+pa[1]-a+d,co[3]+1+b+pa[3]-c,co[4]+1+pa[2]-b+pa[4]-d)
 logw <- rowSums(lgamma(alpha))
 for(j in 1:4) logw <- logw+lchoose(pa[j],z[,j])
 w <- exp(logw-max(logw));w <- w/sum(w);A <- sum(alpha[1,]);H <- co[4]+pa[2]-b+pa[4]-d;mu <- alpha[,4]/A;means <- H+both*mu;m <- sum(w*means)
 list(alloc=z,alpha=alpha,H=H,weights=w,A=A,mean=m,variance=sum(w*(both*mu*(1-mu)*(A+both)/(A+1)+(means-m)^2)),Ep=colSums(alpha*w/A),partialhigh=c(sum(w*(pa[2]-b)),sum(w*(pa[4]-d))))
}
choose <- function(u,p) {
 if(length(u)!=1 || !is.finite(u) || u<=0 || u>=1 || any(!is.finite(p)) || any(p<0) || abs(sum(p)-1)>1e-10) stop('categorical support')
 cum <- cumsum(p)
 if(any(abs(cum[-length(cum)]-u)<=margin)) stop('ambiguous category boundary')
 k <- which(u<cum)[1]
 if(is.na(k)) stop('cumulative failure')
 k
}
p <- read('data/anchor_participants.csv');br <- read('data/anchor_baseline.csv');p <- p[order(p$student_id),]
if(!all(c('student_id','college_id','part_time','first_gen','responded_baseline') %in% names(p)) || !all(c('student_id','college_id','st1','st2') %in% names(br))) stop('required columns')
if(nrow(p)!=1600 || anyDuplicated(p$student_id) || anyDuplicated(br$student_id) || !setequal(br$student_id,p$student_id[p$responded_baseline=='1'])) stop('ID/response inventory')
for(k in c('part_time','first_gen','responded_baseline')) if(any(!p[[k]] %in% c('0','1'))) stop('frame flag')
idx <- match(p$student_id,br$student_id);present <- !is.na(idx)
if(any(br$college_id[idx[present]]!=p$college_id[present])) stop('college join')
raw <- matrix('',1600,2);raw[present,1] <- br$st1[idx[present]];raw[present,2] <- br$st2[idx[present]]
if(any(!raw %in% c('','1','2','3','4','5'))) stop('raw support')
binary <- matrix(-1L,1600,2);binary[raw!=''] <- as.integer(as.numeric(raw[raw!=''])>=4);av <- binary[,1];bv <- binary[,2];group <- paste0(p$part_time,p$first_gen);role <- rep('',1600);role[av<0&bv<0] <- 'both'
for(k in 0:1) {role[av==k&bv<0] <- paste0('row',k);role[av<0&bv==k] <- paste0('col',k)}
fits <- list();groups <- list();totalmean <- totalvar <- known11 <- possible <- bothknown <- 0
for(g in labels) {
 ix <- group==g;if(!any(ix)) stop('empty group')
 co <- sapply(0:3,function(k) sum(ix&av==k%/%2&bv==k%%2));pa <- unname(sapply(rolelabels[1:4],function(k) sum(ix&role==k)));both <- sum(ix&role=='both');f <- mixture(as.integer(co),as.integer(pa),as.integer(both));fits[[g]] <- f;groups[[g]] <- list(counts=co,partial=pa,both=both,roles=setNames(lapply(rolelabels,function(k) which(ix&role==k)),rolelabels));pre <- paste0('group/',g)
 values <- c(n=sum(ix),both_known=sum(co),both_missing=both,configurations=nrow(f$alloc),posterior_count_mean=f$mean,posterior_count_variance=f$variance)
 for(k in names(values)) put(paste0(pre,'/',k),values[[k]])
 for(j in 1:4) {put(paste0(pre,'/observed/',states[j]),co[j]);put(paste0(pre,'/Ep/',states[j]),f$Ep[j]);put(paste0(pre,'/',rolelabels[j]),pa[j])}
 for(j in seq_len(nrow(f$alloc))) put(paste0(pre,'/allocation/',paste0(f$alloc[j,],collapse='_'),'/weight'),f$weights[j])
 put(paste0(pre,'/odds_ratio_at_Ep'),f$Ep[1]*f$Ep[4]/(f$Ep[2]*f$Ep[3]));totalmean <- totalmean+f$mean;totalvar <- totalvar+f$variance;known11 <- known11+co[4];possible <- possible+co[4]+pa[2]+pa[4]+both;bothknown <- bothknown+sum(co)
}
prob <- numeric(1600);eventknown <- av==0|bv==0|(av==1&bv==1)
for(i in 1:1600) {a <- av[i];b <- bv[i];f <- fits[[group[i]]];z <- groups[[group[i]]];prob[i] <- if(a>=0&&b>=0) a*b else if(a==0||b==0) 0 else if(a==1) f$partialhigh[1]/z$partial[2] else if(b==1) f$partialhigh[2]/z$partial[4] else f$Ep[4]}
people <- data.frame(student_id=p$student_id,group=group,raw_st1=raw[,1],raw_st2=raw[,2],observed_A=ifelse(av<0,'',as.character(av)),observed_B=ifelse(bv<0,'',as.character(bv)),event_known=as.character(as.integer(eventknown)),observed_event=ifelse(eventknown,as.character(as.integer(prob)),''),posterior_event_probability=fmt(prob),stringsAsFactors=FALSE)
values <- c('roster/n'=1600,'observed/both_known'=bothknown,'observed/known11'=known11,'support/lower'=known11/1600,'support/upper'=possible/1600,'posterior/mean'=totalmean/1600,'posterior/variance'=totalvar/1600^2,'posterior/SD'=sqrt(totalvar)/1600,'posterior/chebyshev_lower'=max(known11/1600,(totalmean-sqrt(totalvar/.05))/1600),'posterior/chebyshev_upper'=min(possible/1600,(totalmean+sqrt(totalvar/.05))/1600))
for(k in names(values)) put(k,values[[k]])
toys <- read('model-specs/unit06_conditionals.csv');if(!identical(toys$case,c('compatible','incompatible','independent'))) stop('toy inventory')
for(rr in seq_len(nrow(toys))) {
 r <- toys[rr,];a <- as.numeric(r[,c('a_given_b0','a_given_b1')]);b <- as.numeric(r[,c('b_given_a0','b_given_a1')]);kernel(a,b,'AB');oa <- a[2]/(1-a[2])/(a[1]/(1-a[1]));ob <- b[2]/(1-b[2])/(b[1]/(1-b[1]));pre <- paste0('toy/',r$case);laws <- list();put(paste0(pre,'/odds_A'),oa);put(paste0(pre,'/odds_B'),ob);put(paste0(pre,'/compatible'),as.integer(abs(log(oa/ob))<1e-12))
 for(order in c('AB','BA')) {t <- kernel(a,b,order);ss <- stationary(t);laws[[order]] <- ss;for(i in 1:4) {put(paste0(pre,'/',order,'/stationary/',i-1),ss[i]);for(j in 1:4) put(paste0(pre,'/',order,'/transition/',i-1,'/',j-1),t[i,j])};put(paste0(pre,'/',order,'/residual'),max(abs(as.vector(ss%*%t)-ss)))}
 put(paste0(pre,'/order_TV'),sum(abs(laws$AB-laws$BA))/2)
}
one_role <- character();one_key <- character()
for(g in labels) for(k in c('allocation',paste0('gamma',states))) {one_role <- c(one_role,'posterior');one_key <- c(one_key,paste0(g,'/',k))}
for(g in labels) for(r in rolelabels) for(i in groups[[g]]$roles[[r]]) {one_role <- c(one_role,r);one_key <- c(one_key,paste0(g,'/',p$student_id[i]))}
tape <- read('random-inputs/unit06_uniforms.csv');wanted <- unlist(lapply(1:20,function(d) paste(d,one_role,one_key,sep='|')));actual <- paste(tape$imputation,tape$role,tape$key,sep='|');uv <- as.numeric(tape$value)
if(!identical(actual,wanted) || any(!is.finite(uv)) || any(uv<=0|uv>=1)) stop('primitive tape contract')
u <- setNames(uv,actual);getu <- function(d,r,k) unname(u[[paste(d,r,k,sep='|')]])
qs <- numeric();completedrows <- list();worst <- 0
for(d in 1:20) {
 completed <- binary
 for(g in labels) {
  f <- fits[[g]];z <- groups[[g]];ix <- choose(getu(d,'posterior',paste0(g,'/allocation')),f$weights);alloc <- f$alloc[ix,];alpha <- f$alpha[ix,];put(paste0('draw/',d,'/',g,'/allocation_index'),ix-1)
  for(j in 1:4) put(paste0('draw/',d,'/',g,'/allocation/',j-1),alloc[j])
  gamma <- qgamma(sapply(states,function(k) getu(d,'posterior',paste0(g,'/gamma',k))),shape=alpha,scale=1);theta <- gamma/sum(gamma)
  for(j in 1:4) {put(paste0('draw/',d,'/',g,'/gamma/',states[j]),gamma[j]);put(paste0('draw/',d,'/',g,'/theta/',states[j]),theta[j])}
  audit(paste0('chain/',d,'/',g),theta)
  for(order in c('AB','BA')) worst <- max(worst,out[[paste0('chain/',d,'/',g,'/',order,'/maxTV/50')]])
  for(j in 1:4) {
   r <- rolelabels[j];indexes <- z$roles[[r]];left <- alloc[j]
   for(h in seq_along(indexes)) {i <- indexes[h];pr <- left/(length(indexes)-h+1);v <- getu(d,r,paste0(g,'/',p$student_id[i]));if(pr>0&&pr<1&&abs(v-pr)<=margin) stop('ambiguous subset');zero <- as.integer(v<pr);left <- left-zero;completed[i,if(substr(r,1,3)=='row') 2 else 1] <- 1L-zero}
   if(left!=0) stop('subset count failed')
  }
  for(i in z$roles$both) {k <- choose(getu(d,'both',paste0(g,'/',p$student_id[i])),theta)-1;completed[i,] <- c(k%/%2,k%%2)}
 }
 qs[d] <- mean(completed[,1]*completed[,2]);put(paste0('draw/',d,'/Q'),qs[d]);completedrows[[d]] <- data.frame(imputation=as.character(d),student_id=p$student_id,group=group,A=as.character(completed[,1]),B=as.character(completed[,2]),event=as.character(completed[,1]*completed[,2]),stringsAsFactors=FALSE)
}
values <- c('chain/maxTV50'=worst,'chain/mixing_warning'=as.integer(worst>1e-8),'MC/Qbar'=mean(qs),'MC/B'=var(qs),'MC/exact_SE'=sqrt(out[['posterior/variance']]/20),'MC/estimated_SE'=sqrt(var(qs)/20))
for(k in names(values)) put(k,values[[k]])
ex <- which(rowSums(raw!='')==1 & (av==1|bv==1))[1];put('example/posterior_event_probability',prob[ex]);put('example/first_A',as.numeric(completedrows[[1]]$A[ex]));put('example/first_B',as.numeric(completedrows[[1]]$B[ex]))
write.csv(data.frame(key=sort(names(out)),value=fmt(out[sort(names(out))])),output,row.names=FALSE,quote=FALSE)
write.csv(people,sub('\\.csv$','_people.csv',output),row.names=FALSE,quote=FALSE,na='');write.csv(do.call(rbind,completedrows),sub('\\.csv$','_completed.csv',output),row.names=FALSE,quote=FALSE,na='')
cat('UNIT06_BASE_R_OK exact_mixture=true fixed_parameter_kernel=true\n')
