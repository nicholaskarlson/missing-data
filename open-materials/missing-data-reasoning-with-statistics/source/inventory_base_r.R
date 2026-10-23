args<-commandArgs(trailingOnly=TRUE);root<-args[1];out<-args[2]
read<-function(n) read.csv(file.path(root,"data",n),colClasses="character",na.strings=NULL)
p<-read("anchor_participants.csv");b<-read("anchor_baseline.csv");f<-read("anchor_followup.csv");g<-read("anchor_registry.csv");a<-read("anchor_observed_registry.csv")
items<-c(paste0("wb",1:5),paste0("st",1:5));miss<-b[,items]==""
r<-c(invited=nrow(p),responded=sum(p$responded_baseline=="1"),nonresponded=sum(p$responded_baseline=="0"),colleges=length(unique(p$college_id)),baseline_blank_cells=sum(miss),baseline_people_any_blank=sum(rowSums(miss)>0),baseline_stress_complete_people=sum(rowSums(miss[,6:10])==0))
for(arm in c("0","1"))r[paste0("randomized_arm_",arm)]<-sum(b$offer==arm)
for(t in c("3","6","12")){
 w<-f[f$wave_month==t,];r[paste0("wave_",t,"_rows")]<-nrow(w);r[paste0("wave_",t,"_responded")]<-sum(w$responded=="1")
}
for(form in c("1","2","3"))r[paste0("wave_6_form_",form)]<-sum(f$wave_month=="6"&f$responded=="1"&f$form==form)
groups<-list(nonrespondent=a$responded_baseline=="0",control=a$responded_baseline=="1"&a$offer=="0",offered=a$responded_baseline=="1"&a$offer=="1")
for(name in names(groups)){
 d<-a[groups[[name]],];r[paste0(name,"_rows")]<-nrow(d);r[paste0(name,"_known")]<-sum(d$outcome_known=="1");r[paste0(name,"_withdrawals")]<-sum(d$withdrawal_label=="1")
}
r["grade_blank"]<-sum(g$gpa12=="");r["grade_blank_consent"]<-sum(g$gpa12==""&g$exit_type=="consent_withdrawn");r["grade_blank_early_nonconsent"]<-sum(g$gpa12==""&g$exit_type!="consent_withdrawn"&as.numeric(g$exit_month)<4)
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
write.table(data.frame(key=paste0("inventory/",names(r)),value=as.integer(r)),out,sep=",",row.names=FALSE,col.names=TRUE,quote=FALSE)
cat("INVENTORY_BASE_R_OK\n")
