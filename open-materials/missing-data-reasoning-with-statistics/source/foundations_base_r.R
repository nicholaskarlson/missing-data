# Independent base-R observation audit and exhaustive hypothetical mask route.
args <- commandArgs(trailingOnly=TRUE); root <- args[1]; output <- args[2]
read <- function(n) read.csv(file.path(root,n),colClasses='character',check.names=FALSE,na.strings=NULL)
rs <- read('model-specs/foundations_r.csv'); expected <- c(case_id='B8-U02-04-FOUNDATIONS-001',n='1600',groups='00|01|10|11',reverse_item='st4',ordinary_target='finite_roster_mean',toy_n='8',mcar_keep='4',x_keep_0='3',x_keep_1='1',fill_keep_per_group='2')
if(anyDuplicated(rs$name) || !identical(setNames(rs$value,rs$name),expected)) stop('fixed R contract')
p <- read('data/anchor_participants.csv'); b <- read('data/anchor_baseline.csv'); p <- p[order(p$student_id),]
if(nrow(p)!=1600 || anyDuplicated(p$student_id) || anyDuplicated(b$student_id) || !setequal(b$student_id,p$student_id[p$responded_baseline=='1'])) stop('ID/response inventory')
for(k in c('responded_baseline','part_time','first_gen')) if(any(!p[[k]] %in% c('0','1'))) stop('frame flag')
idx <- match(p$student_id,b$student_id); present <- !is.na(idx)
if(!all(c('student_id','college_id',paste0('st',1:5)) %in% names(b)) || !all(c('student_id','college_id','responded_baseline','part_time','first_gen') %in% names(p))) stop('required columns')
for(k in c('college_id')) if(any(b[[k]][idx[present]]!=p[[k]][present])) stop('frame join')
raw <- matrix('',1600,5)
for(j in 1:5) raw[present,j] <- b[[paste0('st',j)]][idx[present]]
if(any(!raw %in% c('','1','2','3','4','5'))) stop('item support')
a <- matrix(NA_real_,1600,5); a[raw!=''] <- as.numeric(raw[raw!='']); a[,4] <- 6-a[,4]
known <- colSums(!is.na(a)); if(any(known==0)) stop('empty item')
itemmeans <- colSums(a,na.rm=TRUE)/known; filled <- a
for(j in 1:5) filled[is.na(a[,j]),j] <- itemmeans[j]
complete <- rowSums(!is.na(a))==5; fs <- rowMeans(filled); groups <- paste0(p$part_time,p$first_gen); out <- numeric()
put <- function(k,v) {if(length(v)!=1 || !is.finite(v)) stop('invalid number');out[k] <<- v}
put('observed/n',1600);put('observed/responded',nrow(b));put('observed/complete',sum(complete));put('observed/complete_mean',mean(rowMeans(a[complete,,drop=FALSE])));put('fill/roster_mean',mean(fs));put('fill/pseudo_iid_variance',var(fs)/1600)
standardized <- 0
for(g in c('00','01','10','11')) {
 z <- groups==g; n <- sum(z); c <- sum(z&complete)
 if(n==0 || c==0) stop('empty complete group')
 m <- mean(rowMeans(a[z&complete,,drop=FALSE])); standardized <- standardized+n/1600*m
 values <- c(n=n,responded=sum(p$responded_baseline[z]=='1'),complete=c,response_rate=mean(p$responded_baseline[z]=='1'),complete_rate=c/n,complete_mean=m,frame_share=n/1600,complete_share=c/sum(complete))
 for(k in names(values)) put(paste0('group/',g,'/',k),values[[k]])
}
put('observed/standardized_complete_mean',standardized)
for(j in 1:5) {put(paste0('item/',j,'/known'),known[j]);put(paste0('item/',j,'/keyed_mean'),itemmeans[j])}
fmt <- function(v) sprintf('%.17g',v)
people <- data.frame(student_id=p$student_id,group=groups,responded=p$responded_baseline,complete_score=ifelse(complete,fmt(rowMeans(a)),''),stringsAsFactors=FALSE)
for(j in 1:5) people[[paste0('raw_st',j)]] <- raw[,j]
for(j in 1:5) people[[paste0('filled_keyed_st',j)]] <- fmt(filled[,j])
people$filled_score <- fmt(fs)
laws <- read('model-specs/foundations_laws.csv');keys <- paste(laws$law,laws$r,laws$y)
if(nrow(laws)!=8 || anyDuplicated(keys) || !setequal(keys,as.vector(outer(c('A','B'),c('0 1','0 5','1 1','1 5'),paste)))) stop('law inventory')
for(label in c('A','B')) {
 z <- laws[laws$law==label,];pr <- as.numeric(z$probability);r <- as.numeric(z$r); y <- as.numeric(z$y);o <- r==1
 if(any(!is.finite(pr)) || any(pr<=0) || abs(sum(pr)-1)>1e-14) stop('law probabilities')
 values <- c(response=sum(pr[o]),observed_mean=sum(y[o]*pr[o])/sum(pr[o]),missing_mean=sum(y[!o]*pr[!o])/sum(pr[!o]),full_mean=sum(y*pr),p_r1_y1=sum(pr[o&y==1])/sum(pr[y==1]),p_r1_y5=sum(pr[o&y==5])/sum(pr[y==5]))
 for(k in names(values)) put(paste0('nonid/',label,'/',k),values[[k]])
 for(i in seq_len(nrow(z))) put(paste0('nonid/',label,'/joint/',r[i],'/',y[i]),pr[i])
}
t <- read('model-specs/foundations_population.csv'); if(nrow(t)!=8 || !identical(t$id,paste0('H',1:8)) || !identical(t$x,as.character(rep(0:1,each=4))) || !identical(t$y,as.character(c(1,1,3,3,3,3,5,5)))) stop('fixed teaching roster')
y <- as.numeric(t$y);x <- as.numeric(t$x);put('toy/true_mean',mean(y));put('toy/true_slope',mean(y[x==1])-mean(y[x==0]))
cross <- function(k0,k1) {
 a <- combn(1:4,k0,simplify=FALSE);b <- combn(5:8,k1,simplify=FALSE);result <- list()
 for(i in a) for(j in b) result[[length(result)+1]] <- c(i,j)
 result
}
masks <- list(mcar=combn(1:8,4,simplify=FALSE),xmask=cross(3,1),fill=cross(2,2))
for(label in names(masks)) {
 measures <- list()
 for(ix in masks[[label]]) {
  v <- y[ix];m <- mean(v);values <- c(mean=m);tag <- paste0(as.integer((1:8) %in% ix),collapse='')
  if(label=='xmask') {m0 <- mean(v[x[ix]==0]);m1 <- mean(v[x[ix]==1]);values <- c(values,standardized=(m0+m1)/2,slope=m1-m0)}
  if(label=='fill') {z <- rep(m,8);z[ix] <- v;values <- c(values,filled_mean=mean(z),retained_iid_variance=var(v)/4,pseudo_iid_variance=var(z)/8,variance_factor=(length(v)-1)/(length(y)-1)*length(v)/length(y))}
  for(k in names(values)) {put(paste0('toy/',label,'/mask/',tag,'/',k),values[[k]]);measures[[k]] <- c(measures[[k]],values[[k]])}
 }
 put(paste0('toy/',label,'/masks'),length(masks[[label]]))
 for(k in names(measures)) {v <- measures[[k]];put(paste0('toy/',label,'/expectation/',k),mean(v));put(paste0('toy/',label,'/mask_variance/',k),mean((v-mean(v))^2))}
}
write.csv(data.frame(key=sort(names(out)),value=fmt(out[sort(names(out))])),output,row.names=FALSE,quote=FALSE)
write.csv(people,sub('\\.csv$','_people.csv',output),row.names=FALSE,quote=FALSE,na='')
cat('FOUNDATIONS_BASE_R_OK\n')
