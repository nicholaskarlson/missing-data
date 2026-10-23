# Independent base-R joins, item arithmetic, support and randomization variance bound.
args <- commandArgs(trailingOnly=TRUE); root <- args[1]; output <- args[2]
read <- function(n) read.csv(file.path(root,n),colClasses='character',na.strings=NULL,check.names=FALSE)
spec <- read('model-specs/unit09_r.csv'); cfg <- setNames(spec$value,spec$name); n_expected <- as.integer(cfg['n']); alpha <- as.numeric(cfg['alpha']); items <- paste0('wb',1:5)
stopifnot(is.finite(alpha),alpha>0,alpha<1,cfg['reverse_item']=='wb4',cfg['target']=='finite_randomized_roster_offer_effect')
p <- read('data/anchor_participants.csv'); b <- read('data/anchor_baseline.csv'); f <- read('data/anchor_followup.csv')
stopifnot(!anyDuplicated(p$student_id),all(p$responded_baseline %in% c('0','1')),all(p$part_time %in% c('0','1')),nrow(b)==n_expected,!anyDuplicated(b$student_id),setequal(b$student_id,p$student_id[p$responded_baseline=='1']))
b <- b[order(b$student_id),]; p <- p[match(b$student_id,p$student_id),]; stopifnot(all(b$college_id==p$college_id),all(b$offer %in% c('0','1')))
keys <- paste(f$student_id,f$wave_month,sep='/'); wanted <- as.vector(outer(b$student_id,c('3','6','12'),paste,sep='/')); stopifnot(!anyDuplicated(keys),setequal(keys,wanted),all(f$responded %in% c('0','1')),all(as.matrix(f[items]) %in% c('',as.character(1:5))),all(as.matrix(f[f$responded=='0',items])==''))
waves <- lapply(c('3','6','12'),function(t) f[match(paste(b$student_id,t,sep='/'),keys),]); last <- waves[[3]]; raw <- as.matrix(last[items]); present <- raw!=''; x <- matrix(0,nrow(b),5); x[present] <- as.numeric(raw[present]); x[,4] <- ifelse(present[,4],6-x[,4],0)
k <- unname(rowSums(present)); total <- rowSums(x); missing <- 5-k; lower <- (total+missing)/5; upper <- (total+5*missing)/5; full <- k==5; score <- total/5; z <- as.integer(b$offer); pt <- as.integer(p$part_time); pattern <- paste0(waves[[1]]$responded,waves[[2]]$responded,waves[[3]]$responded); n <- nrow(b); colleges <- sort(unique(b$college_id)); stopifnot(length(colleges)==as.integer(cfg['college_count']))
out <- c(randomized=n,college_count=length(colleges))
group <- function(m) {
 nn <- sum(m); stopifnot(nn>0); q <- m&full; nc <- sum(q)
 a <- c(n=nn,responded12=sum(last$responded[m]=='1'),complete=nc,known_cells=sum(k[m]),keyed_sum=sum(total[m]),part_time=sum(pt[m]),complete_part_time=sum(pt[q]),lower=mean(lower[m]),upper=mean(upper[m]),complete_case_defined=as.integer(nc>0))
 a <- c(a,missing_cells=5*nn-a['known_cells'][[1]],response_fraction=a['responded12'][[1]]/nn,complete_fraction=nc/nn,part_time_fraction=a['part_time'][[1]]/nn)
 if(nc>0) a <- c(a,complete_mean=mean(score[q]),complete_part_time_fraction=mean(pt[q]))
 a
}
put <- function(prefix,a) {names(a)<-paste0(prefix,names(a));out <<- c(out,a)}
put('all/',group(rep(TRUE,n))); arms <- lapply(0:1,function(a) group(z==a)); for(a in 0:1) put(paste0('arm/',a,'/'),arms[[a+1]])
ok <- all(vapply(arms,function(g)g['complete_case_defined']==1,logical(1))); out['complete_case/pooled_defined'] <- as.integer(ok)
if(ok) out['complete_case/pooled_difference'] <- unname(arms[[2]]['complete_mean']-arms[[1]]['complete_mean'])
for(pa in c('000','001','010','011','100','101','110','111')) for(a in 0:1) out[paste0('pattern/',pa,'/',a)] <- sum(pattern==pa&z==a)
out['intermittent_return'] <- sum(grepl('01',pattern,fixed=TRUE)); dl<-du<-vb<-cc<-numeric(); defined<-TRUE
for(college in colleges) {
 mask<-b$college_id==college; nc<-sum(mask); w<-nc/n; g0<-group(mask&z==0); g1<-group(mask&z==1); stopifnot(nc>=2,abs(g1['n']-g0['n'])<=1)
 prefix<-paste0('college/',college,'/'); put(prefix,c(n=nc,weight=w));put(paste0(prefix,'0/'),g0);put(paste0(prefix,'1/'),g1)
 low<-unname(g1['lower']-g0['upper']);high<-unname(g1['upper']-g0['lower']); v<-unname(4*w*w*nc/(nc-1)*(1/g1['n']+1/g0['n']));put(prefix,c(contrast_lower=low,contrast_upper=high,variance_bound=v));dl<-c(dl,w*low);du<-c(du,w*high);vb<-c(vb,v)
 ok<-g0['complete_case_defined']==1&&g1['complete_case_defined']==1;defined<-defined&&ok;if(ok)cc<-c(cc,w*(g1['complete_mean']-g0['complete_mean']))
}
out['complete_case/standardized_defined']<-as.integer(defined);if(defined)out['complete_case/standardized_difference']<-sum(cc)
low<-sum(dl);high<-sum(du);vmax<-sum(vb);radius<-sqrt(vmax/alpha); put('realized/',c(lower=low,upper=high,width=high-low));put('confidence/',c(alpha=alpha,variance_bound=vmax,radius=radius,lower=max(-4,low-radius),upper=min(4,high+radius)));put('causal_support/',c(lower=mean(ifelse(z==1,lower-5,1-upper)),upper=mean(ifelse(z==1,upper-1,5-lower))))
example<-which(last$responded=='1'&missing==1&present[,4]&raw[,4]!='3');out['example/exists']<-as.integer(length(example)>0)
if(length(example)>0) {i<-example[1];put('example/',c(observed_items=k[i],keyed_sum=total[i],lower=lower[i],upper=upper[i],missing_item=unname(which(!present[i,])),keyed_wb4=x[i,4]));for(j in 1:5)if(present[i,j])out[paste0('example/raw/',items[j])]<-as.numeric(raw[i,j])}
stopifnot(all(is.finite(out)),!anyDuplicated(names(out)));out<-out[order(names(out))];write.table(data.frame(key=names(out),value=sprintf('%.17g',out)),output,sep=',',quote=FALSE,row.names=FALSE)
people<-data.frame(student_id=b$student_id,college_id=b$college_id,offer=b$offer,part_time=p$part_time,wave_pattern=pattern,responded12=last$responded,observed_items=k,keyed_sum=total,complete_score=ifelse(full,sprintf('%.17g',score),''),lower=sprintf('%.17g',lower),upper=sprintf('%.17g',upper));write.table(people,sub('\\.csv$','_people.csv',output),sep=',',quote=FALSE,row.names=FALSE)
cat('UNIT09_BASE_R_OK randomization_bound=true response_model=none\n')
