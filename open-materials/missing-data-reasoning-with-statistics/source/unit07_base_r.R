# Independent base-R pattern likelihood, EM and information calculation.
args<-commandArgs(TRUE);root<-args[1];outfile<-args[2]
read<-function(p)read.csv(file.path(root,p),colClasses='character',na.strings=NULL,check.names=FALSE)
spec<-read('model-specs/unit07_r.csv');stopifnot(!anyDuplicated(spec$name));s<-setNames(spec$value,spec$name)
expected<-c(case_id='B8-U07-OBSERVED-LIKELIHOOD',n='1600',em_sweeps='100',prior='none',groups='00|01|10|11',high_threshold='4',scaled_score_limit='1e-12',step_limit='1e-13',likelihood_decrease_limit='1e-10',information_floor='1e-12',newton_sweeps='20',wald_z='1.959963984540054')
stopifnot(identical(sort(names(s)),sort(names(expected))),all(s[names(expected)]==expected))
frame<-read('data/anchor_participants.csv');base<-read('data/anchor_baseline.csv');frame<-frame[order(frame$student_id),];stopifnot(nrow(frame)==1600,!anyDuplicated(frame$student_id),!anyDuplicated(base$student_id),setequal(base$student_id,frame$student_id[frame$responded_baseline=='1']))
stopifnot(all(frame$part_time%in%c('0','1')),all(frame$first_gen%in%c('0','1')),all(frame$responded_baseline%in%c('0','1')))
join<-match(frame$student_id,base$student_id);has<-!is.na(join);stopifnot(all(frame$college_id[has]==base$college_id[join[has]]))
raw<-matrix('',1600,2);raw[has,]<-as.matrix(base[join[has],c('st1','st2')]);stopifnot(all(raw%in%c('','1','2','3','4','5')))
bits<-matrix(NA_integer_,1600,2);present<-raw!='';bits[present]<-as.integer(raw[present])>=4
states<-c('00','01','10','11');roles<-c('row0','row1','col0','col1');group<-paste0(frame$part_time,frame$first_gen)
mask<-rbind(diag(4),c(1,1,0,0),c(0,0,1,1),c(1,0,1,0),c(0,1,0,1));design<-mask[,2:4]-mask[,1]
pattern<-character(1600);category<-integer(1600)
for(i in 1:1600){a<-bits[i,1];b<-bits[i,2];if(is.na(a)&&is.na(b)){pattern[i]<-'neither';category[i]<-9}else if(is.na(a)){pattern[i]<-paste0('col',b);category[i]<-7+b}else if(is.na(b)){pattern[i]<-paste0('row',a);category[i]<-5+a}else{pattern[i]<-'complete';category[i]<-1+2*a+b}}
prob<-function(p){stopifnot(length(p)==4,all(is.finite(p)),all(p>0),abs(sum(p)-1)<1e-12);as.vector(mask%*%p)}
ll<-function(cnt,p)sum(cnt*log(prob(p)))
curv<-function(cnt,p){q<-prob(p);list(score=as.vector(crossprod(design,cnt/q)),info=crossprod(design,design*(cnt/q^2)))}
inv<-function(h){e<-eigen(h,symmetric=TRUE,only.values=TRUE)$values;stopifnot(all(is.finite(h)),min(e)>1e-12*max(e));solve(h)}
em<-function(cnt,start,principal=TRUE){stopifnot(length(cnt)==8,all(cnt>=0),all(cnt==as.integer(cnt)),sum(cnt)>0);if(principal)stopifnot(all(cnt[1:4]>0));p<-start;value<-ll(cnt,p);trace<-list('0'=list(p=p,ll=value));decrease<-0
 for(it in 1:100){allocation<-sweep(mask,2,p,'*');allocation<-allocation/as.vector(mask%*%p);pn<-as.vector(crossprod(allocation,cnt))/sum(cnt);now<-ll(cnt,pn);decrease<-max(decrease,value-now);stopifnot(now>=value-1e-10);last<-max(abs(pn-p));p<-pn;value<-now;if(it%in%c(1,2,5,10,25,50,100))trace[[as.character(it)]]<-list(p=p,ll=value)}
 if(principal)stopifnot(last<=1e-13,max(abs(curv(cnt,p)$score))/sum(cnt)<=1e-12)
 list(p=p,trace=trace,last=last,decrease=decrease)}
newton<-function(cnt){p<-rep(.25,4);for(it in 1:20){d<-curv(cnt,p);step<-solve(d$info,d$score);old<-ll(cnt,p);found<-FALSE;for(h in 0:40){t<-p[2:4]+step*2^(-h);nextp<-c(1-sum(t),t);if(min(nextp)>0&&ll(cnt,nextp)>=old-1e-10){found<-TRUE;break}};stopifnot(found);p<-nextp};stopifnot(max(abs(curv(cnt,p)$score))/sum(cnt)<=1e-12);p}
output<-list();put<-function(k,v){stopifnot(length(v)==1,is.finite(v));output[[k]]<<-as.numeric(v)}
fits<-list();estimate<-variance<-naive<-cc<-maxspread<-0
starts<-list(rep(.25,4),c(.7,.1,.1,.1),c(.1,.2,.3,.4))
for(g in states){sel<-group==g;counts<-tabulate(category[sel],9);cnt<-counts[1:8];stopifnot(all(cnt[1:4]>0));prefix<-paste0('group/',g);solutions<-list()
 for(st in 0:2){fit<-em(cnt,starts[[st+1]]);solutions[[st+1]]<-fit$p;put(paste0(prefix,'/start/',st,'/last_step'),fit$last);put(paste0(prefix,'/start/',st,'/max_likelihood_decrease'),fit$decrease)
  for(it in names(fit$trace)){tr<-fit$trace[[it]];put(paste0(prefix,'/start/',st,'/sweep/',it,'/loglik'),tr$ll);for(j in 1:4)put(paste0(prefix,'/start/',st,'/sweep/',it,'/p',states[j]),tr$p[j])}}
 p<-solutions[[1]];fits[[g]]<-p;d<-curv(cnt,p);cov<-inv(d$info);direct<-newton(cnt);spread<-max(abs(unlist(solutions)-rep(p,3)));maxspread<-max(maxspread,spread);stopifnot(max(abs(p-direct))<=1e-11)
 put(paste0(prefix,'/n'),sum(counts));put(paste0(prefix,'/informative'),sum(cnt));put(paste0(prefix,'/neither'),counts[9]);put(paste0(prefix,'/loglik'),ll(cnt,p));put(paste0(prefix,'/scaled_score'),max(abs(d$score))/sum(cnt));put(paste0(prefix,'/start_spread'),spread)
 for(j in 1:4){put(paste0(prefix,'/complete/',states[j]),cnt[j]);put(paste0(prefix,'/p',states[j]),p[j]);put(paste0(prefix,'/newton/p',states[j]),direct[j]);put(paste0(prefix,'/',roles[j]),cnt[4+j])}
 for(i in 1:3){put(paste0(prefix,'/score/',i-1),d$score[i]);for(j in 1:3){put(paste0(prefix,'/information/',i-1,'/',j-1),d$info[i,j]);put(paste0(prefix,'/covariance/',i-1,'/',j-1),cov[i,j])}}
 weight<-sum(counts)/1600;estimate<-estimate+weight*p[4];variance<-variance+weight^2*cov[3,3];naive<-naive+weight^2*p[4]*(1-p[4])/sum(counts);cc<-cc+weight*cnt[4]/sum(cnt[1:4]);put(paste0(prefix,'/observed_SE_p11'),sqrt(cov[3,3]));put(paste0(prefix,'/fractional_all_rows_SE_p11'),sqrt(p[4]*(1-p[4])/sum(counts)))}
rows<-data.frame(student_id=frame$student_id,group=group,raw_st1=raw[,1],raw_st2=raw[,2],observed_A=ifelse(is.na(bits[,1]),'',bits[,1]),observed_B=ifelse(is.na(bits[,2]),'',bits[,2]),pattern=pattern,event_known='',observed_event='',plugin_event_probability='',stringsAsFactors=FALSE)
probabilities<-numeric(1600);lower<-upper<-0
for(i in 1:1600){a<-bits[i,1];b<-bits[i,2];p<-fits[[group[i]]];known<-(!is.na(a)&&a==0)||(!is.na(b)&&b==0)||(!is.na(a)&&!is.na(b));v<-if(!is.na(a)&&!is.na(b))a*b else if((!is.na(a)&&a==0)||(!is.na(b)&&b==0))0 else if(!is.na(a))p[4]/(p[3]+p[4]) else if(!is.na(b))p[4]/(p[2]+p[4]) else p[4];probabilities[i]<-v;rows$event_known[i]<-as.character(as.integer(known));if(known)rows$observed_event[i]<-as.character(as.integer(v));rows$plugin_event_probability[i]<-sprintf('%.17g',v);lower<-lower+as.integer(!is.na(a)&&!is.na(b)&&a==1&&b==1);upper<-upper+as.integer((is.na(a)||a!=0)&&(is.na(b)||b!=0))}
se<-sqrt(variance);put('model/standardized_probability',estimate);put('model/observed_variance',variance);put('model/observed_SE',se);put('model/asymptotic_wald_lower',estimate-as.numeric(s['wald_z'])*se);put('model/asymptotic_wald_upper',estimate+as.numeric(s['wald_z'])*se);put('naive/fractional_all_rows_SE',sqrt(naive));put('description/complete_pair_standardized',cc);put('finite/plugin_prediction',mean(probabilities));put('finite/logical_lower',lower/1600);put('finite/logical_upper',upper/1600);put('audit/max_start_spread',maxspread);put('audit/prediction_score_identity_residual',abs(mean(probabilities)-estimate))
tables<-list(c(.4,.1,.1,.4),c(.1,.4,.4,.1));cnt<-c(rep(0,4),rep(10,4))
for(k in 0:1){p<-tables[[k+1]];put(paste0('toy/ridge/',k,'/loglik'),ll(cnt,p));fit<-em(cnt,p,FALSE);put(paste0('toy/ridge/',k,'/stationary_p11'),fit$p[4]);h<-curv(cnt,p);put(paste0('toy/ridge/',k,'/scaled_score'),max(abs(h$score))/sum(cnt));fail<-tryCatch({inv(h$info);FALSE},error=function(e)TRUE);stopifnot(fail);put(paste0('toy/ridge/',k,'/information_singular'),1);put(paste0('toy/no_observations/',k,'/loglik'),ll(rep(0,8),p))}
cnt<-c(4,1,1,4,rep(0,4));fit<-em(cnt,rep(.25,4));put('toy/complete/p11',fit$p[4]);put('toy/complete/parameter_variance',inv(curv(cnt,fit$p)$info)[3,3]);put('toy/complete/finite_missing_variance',0)
keys<-sort(names(output));write.table(data.frame(key=keys,value=vapply(output[keys],function(v)sprintf('%.17g',v),character(1))),outfile,sep=',',row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(rows,sub('\\.csv$','_people.csv',outfile),sep=',',row.names=FALSE,col.names=TRUE,quote=FALSE,na='');cat('UNIT07_BASE_R_OK observed_likelihood=true finite_interval=false\n')
