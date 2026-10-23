# Independent base-R row completion and college standardization; no analysis packages.
args<-commandArgs(trailingOnly=TRUE);root<-args[1];output<-args[2]
rd<-function(n)read.csv(file.path(root,n),colClasses='character',check.names=FALSE,na.strings=NULL)
need<-function(x,msg)if(!isTRUE(x))stop(msg)
items<-paste0('wb',1:5);grid<-c(-8,-4,-2,-1,0,1,2,4,8)
cs<-rd('model-specs/unit11_r.csv');need(!anyDuplicated(cs$name),'duplicate contract key');s<-setNames(cs$value,cs$name)
need(identical(as.character(s[c('n','college_count','wave','items','reverse_item','support','grid_half_points','scenario_half_points','paths')]),c('1219','16','12','wb1|wb2|wb3|wb4|wb5','wb4','1|5','-8|-4|-2|-1|0|1|2|4|8','1|-1','offered_worse|differential')),'fixed Unit11 R specification')
f<-rd('data/anchor_participants.csv');b<-rd('data/anchor_baseline.csv');v<-rd('data/anchor_followup.csv')
need(nrow(f)==1600&&!anyDuplicated(f$student_id)&&all(f$responded_baseline%in%c('0','1')),'frame identity/flags')
need(nrow(b)==1219&&!anyDuplicated(b$student_id)&&setequal(b$student_id,f$student_id[f$responded_baseline=='1']),'randomized membership')
need(all(b$offer%in%c('0','1'))&&all(b$college_id==f$college_id[match(b$student_id,f$student_id)]),'offer/college join')
keys<-paste(v$student_id,v$wave_month);wanted<-as.vector(outer(b$student_id,c('3','6','12'),paste));need(!anyDuplicated(keys)&&setequal(keys,wanted),'follow-up inventory')
need(all(v$responded%in%c('0','1'))&&all(as.matrix(v[items])%in%c('',as.character(1:5))),'item/response support')
need(all(as.matrix(v[v$responded=='0',items])==''),'wave nonresponse contains answers')
b<-b[order(b$student_id),];last<-v[v$wave_month=='12',];last<-last[match(b$student_id,last$student_id),]
p<-data.frame(student_id=b$student_id,college_id=b$college_id,offer=b$offer,responded12=last$responded,last[items],stringsAsFactors=FALSE)
need(length(unique(p$college_id))==16,'college count')
for(c in unique(p$college_id)){n<-table(factor(p$offer[p$college_id==c],levels=c('0','1')));need(all(n>0)&&abs(n[1]-n[2])<=1,'original balance')}
toy<-rd('model-specs/unit11_teaching.csv');need(nrow(toy)==12&&setequal(toy$fixture,c('gaps','complete'))&&!anyDuplicated(paste(toy$fixture,toy$id))&&all(toy$offer%in%c('0','1'))&&all(as.matrix(toy[items])%in%c('',as.character(1:5))),'teaching fixture')
prepare<-function(p){
 X<-matrix(as.numeric(as.matrix(p[items])),nrow=nrow(p),ncol=5);X[,4]<-6-X[,4];A<-X;cells<-list();index<-1
 for(c in sort(unique(p$college_id))){nc<-sum(p$college_id==c)
  for(z in c('0','1')){ix<-which(p$college_id==c&p$offer==z);n<-length(ix);need(n>0,'empty arm')
   for(j in 1:5){known<-X[ix,j][!is.na(X[ix,j])];need(length(known)>0,'unavailable item anchor');anchor<-sum(known)/length(known);A[ix,j]<-anchor
    cells[[index]]<-data.frame(college_id=c,offer=z,item=items[j],n=n,college_n=nc,observed=length(known),missing=n-length(known),wave_missing=sum(p$responded12[ix]=='0'),item_missing=sum(is.na(X[ix,j])&p$responded12[ix]=='1'),keyed_sum=sum(known),anchor=anchor,weight=nc/nrow(p),stringsAsFactors=FALSE);index<-index+1
   }
  }
 }
 list(p=p,X=X,A=A,cells=do.call(rbind,cells))
}
estimate<-function(q,d0,d1){
 D<-ifelse(q$p$offer=='0',d0,d1);pred<-pmin(5,pmax(1,q$A+D));full<-q$X;full[is.na(full)]<-pred[is.na(full)];score<-rowMeans(full);m<-c(0,0)
 for(c in sort(unique(q$p$college_id))){w<-sum(q$p$college_id==c)/nrow(q$p)
  for(z in 0:1)m[z+1]<-m[z+1]+w*mean(score[q$p$college_id==c&q$p$offer==as.character(z)])
 }
 c(control=m[1],offered=m[2],D=m[2]-m[1])
}
first_zero<-function(q,path){
 slope<-if(path=='offered_worse')c(0,-1) else c(.5,-.5);limit<-if(path=='offered_worse')4 else 8
 fun<-function(t)unname(estimate(q,slope[1]*t,slope[2]*t)['D']);start<-fun(0);end<-fun(limit);out<-c(start=start,end=end,already_nonpositive=as.numeric(start<=0),defined=0)
 if(start<0||end>0)return(out)
 if(start==0)return(c(out[c('start','end','already_nonpositive')],defined=1,t=0,residual=0))
 lo<-0;hi<-limit
 for(i in 1:80){mid<-(lo+hi)/2;if(fun(mid)>0)lo<-mid else hi<-mid}
 c(out[c('start','end','already_nonpositive')],defined=1,t=hi,residual=fun(hi))
}
summarize<-function(p){
 q<-prepare(p);cells<-q$cells;o<-c(n=nrow(p),cells=nrow(cells),known_cells=sum(cells$observed),missing_cells=sum(cells$missing),wave_missing_cells=sum(cells$wave_missing),item_missing_cells=sum(cells$item_missing))
 for(z in c('0','1'))o<-c(o,setNames(c(sum(p$offer==z),sum(cells$missing[cells$offer==z])),paste0('arm/',z,c('/n','/missing_cells'))))
 scenarios<-list(reference=c(0,0),scenario=c(.5,-.5),support_lower=c(4,-4),support_upper=c(-4,4))
 for(label in names(scenarios)){d<-scenarios[[label]];a<-estimate(q,d[1],d[2]);o<-c(o,setNames(a,paste0(label,'/',names(a))))}
 gg<-list();k<-1
 for(i in seq_along(grid))for(j in seq_along(grid)){
  d0<-grid[i]/2;d1<-grid[j]/2;a<-estimate(q,d0,d1);shift<-ifelse(cells$offer=='0',d0,d1);lower<-sum(cells$missing[cells$anchor+shift<=1]);upper<-sum(cells$missing[cells$anchor+shift>=5]);row<-data.frame(control_half_points=grid[i],offered_half_points=grid[j],control_mean=unname(a[1]),offered_mean=unname(a[2]),D=unname(a[3]),lower_boundary_cells=lower,upper_boundary_cells=upper);gg[[k]]<-row;k<-k+1
  vals<-unlist(row[c('control_mean','offered_mean','D','lower_boundary_cells','upper_boundary_cells')]);o<-c(o,setNames(vals,paste0('grid/',i-1,'/',j-1,'/',names(vals))))
 }
 for(path in c('offered_worse','differential')){a<-first_zero(q,path);o<-c(o,setNames(a,paste0('tip/',path,'/',names(a))))}
 g<-do.call(rbind,gg);common<-g$D[g$control_half_points==g$offered_half_points];o<-c(o,'common/grid_min_D'=min(common),'common/grid_max_D'=max(common))
 fill<-function(d){pred<-pmin(5,pmax(1,q$A+d));x<-q$X;x[is.na(x)]<-pred[is.na(x)];rowMeans(x)}
 people<-data.frame(p,observed_items=rowSums(!is.na(q$X)),keyed_sum=rowSums(q$X,na.rm=TRUE),reference_score=fill(rep(0,nrow(p))),scenario_score=fill(ifelse(p$offer=='0',.5,-.5)))
 list(o=o,people=people,cells=cells,grid=g,q=q)
}
a<-summarize(p);o<-a$o
for(fixture in c('gaps','complete')){
 t<-toy[toy$fixture==fixture,];pp<-data.frame(student_id=t$id,college_id=t$college_id,offer=t$offer,responded12=ifelse(rowSums(t[items]=='')==5,'0','1'),t[items]);v<-summarize(pp)
 for(k in c('reference/D','scenario/D','support_lower/D','support_upper/D'))o<-c(o,setNames(v$o[k],paste0('toy/',fixture,'/',k)))
 for(path in c('offered_worse','differential')){x<-first_zero(v$q,path);o<-c(o,setNames(x,paste0('toy/',fixture,'/tip/',path,'/',names(x))))}
}
need(all(is.finite(o))&&!anyDuplicated(names(o)),'nonfinite/duplicate output');options(digits=17)
write.table(data.frame(key=names(o)[order(names(o))],value=sprintf('%.17g',o[order(names(o))])),output,sep=',',row.names=FALSE,quote=FALSE)
for(name in c('people','cells','grid'))write.table(a[[name]],sub('\\.csv$',paste0('_',name,'.csv'),output),sep=',',row.names=FALSE,quote=FALSE,na='')
cat('UNIT11_BASE_R_OK bounded_mean_sensitivity=true interval=false\n')
