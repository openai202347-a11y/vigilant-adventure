args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2) stop('Usage: Rscript run.R PAPER_DIR R_LIBRARY')
root <- normalizePath(args[1],winslash='/',mustWork=TRUE)
.libPaths(c(normalizePath(args[2],winslash='/',mustWork=TRUE),.libPaths()))
suppressPackageStartupMessages(library(survey)); suppressPackageStartupMessages(library(foreign))
options(survey.lonely.psu='fail',survey.adjust.domain.lonely=FALSE)
raw <- file.path(root,'03_data','raw'); out <- file.path(root,'04_analysis','results'); dir.create(out,recursive=TRUE,showWarnings=FALSE)
readpart <- function(name,cols) { x<-foreign::read.xport(file.path(raw,paste0(name,'.xpt'))); names(x)<-toupper(names(x)); stopifnot(all(cols %in% names(x)),!anyDuplicated(x$SEQN)); x[,cols,drop=FALSE] }
d <- readpart('DEMO_L',c('SEQN','RIDAGEYR','RIAGENDR','RIDEXPRG','RIDRETH3','SDMVSTRA','SDMVPSU'))
parts<-list(DIQ_L=c('SEQN','DIQ010'),GHB_L=c('SEQN','LBXGH','WTPH2YR'),BMX_L=c('SEQN','BMXBMI'),BPXO_L=c('SEQN','BPXOSY1','BPXOSY2','BPXOSY3','BPXODI1','BPXODI2','BPXODI3'),BPQ_L=c('SEQN','BPQ020','BPQ150'))
for(n in names(parts)) d<-merge(d,readpart(n,parts[[n]]),by='SEQN',all.x=TRUE,sort=FALSE)
stopifnot(!anyDuplicated(d$SEQN)); avg<-function(x){z<-rowMeans(x,na.rm=TRUE); z[!is.finite(z)]<-NA_real_; z}
d$sys<-avg(d[,c('BPXOSY1','BPXOSY2','BPXOSY3')]);d$dia<-avg(d[,c('BPXODI1','BPXODI2','BPXODI3')])
d$highbp<-NA_integer_; pos<-(!is.na(d$sys)&d$sys>=140)|(!is.na(d$dia)&d$dia>=90)|(!is.na(d$BPQ150)&d$BPQ150==1)
neg<-(!is.na(d$sys)&d$sys<140)&(!is.na(d$dia)&d$dia<90)&((!is.na(d$BPQ020)&d$BPQ020==2)|(!is.na(d$BPQ150)&d$BPQ150==2))
d$highbp[pos]<-1L; d$highbp[!pos & neg]<-0L
age<-d$RIDAGEYR; bmi<-d$BMXBMI
d$agepoints<-ifelse(age>=25 & age<35,0,ifelse(age<45,8,ifelse(age<55,13,ifelse(age<65,20,ifelse(age<=74,22,NA)))))
d$bmipoints<-ifelse(is.na(bmi)|bmi<18,NA,ifelse(bmi<20,-6,ifelse(bmi<23,0,ifelse(bmi<25,9,ifelse(bmi<27.5,14,ifelse(bmi<30,18,20))))))
d$score<-d$agepoints+d$bmipoints+8*d$highbp; d$case<-as.integer(d$LBXGH>=6.5);d$positive<-as.integer(d$score>=28)
flow<-data.frame(stage=character(),remaining=integer(),excluded=integer())
keep<-rep(TRUE,nrow(d)); step<-function(label,condition){condition[is.na(condition)]<-FALSE; previous<-sum(keep); keep<<-keep & condition; flow<<-rbind(flow,data.frame(stage=label,remaining=sum(keep),excluded=previous-sum(keep)))}
step('All interview records',rep(TRUE,nrow(d)));step('Age 25-74',age>=25&age<=74)
step('No self-reported diabetes (DIQ010=2)',d$DIQ010==2)
pregnant_restricted<-d$RIAGENDR==2 & age>=25 & age<=44
step('Pregnancy: exclude known and unknown among women 25-44',!pregnant_restricted | d$RIDEXPRG==2)
step('HbA1c available',!is.na(d$LBXGH));step('BMI >=18 and available',!is.na(bmi)&bmi>=18);step('NDS blood pressure component available',!is.na(d$highbp))
step('Valid sampling design and blood weight',!is.na(d$WTPH2YR)&d$WTPH2YR>0&!is.na(d$SDMVSTRA)&!is.na(d$SDMVPSU))
d$eligible<-keep
write.csv(flow,file.path(out,'flow.csv'),row.names=FALSE,fileEncoding='UTF-8')
# Survey design retains the complete sampled blood-weight population; domains are defined after design construction.
all<-d[!is.na(d$WTPH2YR)&d$WTPH2YR>0&!is.na(d$SDMVSTRA)&!is.na(d$SDMVPSU),]
des<-svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~WTPH2YR,nest=TRUE,data=all)
main<-subset(des,eligible); x<-main$variables
if(nrow(x)==0 || length(unique(x$case))<2 || length(unique(x$positive))<2) stop('Insufficient cases or score classes; do not report metrics')
metric<-function(design,domain,variable,label){q<-subset(design,eval(parse(text=domain))); if(nrow(q$variables)<2) return(data.frame(metric=label,n=nrow(q$variables),estimate=NA,lower=NA,upper=NA)); a<-tryCatch(svymean(as.formula(paste0('~',variable)),q,na.rm=TRUE),error=function(e) NULL);if(is.null(a)) return(data.frame(metric=label,n=nrow(q$variables),estimate=NA,lower=NA,upper=NA)); ci<-confint(a); data.frame(metric=label,n=nrow(q$variables),estimate=unname(coef(a)[1]),lower=max(0,ci[1,1]),upper=min(1,ci[1,2]))}
metrics<-do.call(rbind,list(metric(main,'case==1','positive','Sensitivity'),metric(main,'case==0','I(1-positive)','Specificity'),metric(main,'positive==1','case','Positive predictive value'),metric(main,'positive==0','I(1-case)','Negative predictive value'),metric(main,'TRUE','positive','Screen-positive proportion'),metric(main,'TRUE','case','Elevated HbA1c prevalence')))
write.csv(metrics,file.path(out,'screening_metrics.csv'),row.names=FALSE)
# Weighted empirical AUC, tie-corrected; no AUC CI if design-consistent resampling is unavailable.
posdat<-x[x$case==1,]; negdat<-x[x$case==0,]; weighted_auc<-sum(vapply(seq_len(nrow(posdat)),function(i) posdat$WTPH2YR[i]*sum(negdat$WTPH2YR*((posdat$score[i]>negdat$score)+0.5*(posdat$score[i]==negdat$score))),numeric(1)))/(sum(posdat$WTPH2YR)*sum(negdat$WTPH2YR))
write.csv(data.frame(metric='Weighted AUC',estimate=weighted_auc,lower=NA,upper=NA,n=nrow(x),events=sum(x$case),note='95% CI not estimated; design-consistent resampling not implemented'),file.path(out,'auc.csv'),row.names=FALSE)
groups<-list('Men'='RIAGENDR==1','Women'='RIAGENDR==2','Age 25-44'='RIDAGEYR>=25 & RIDAGEYR<=44','Age 45-64'='RIDAGEYR>=45 & RIDAGEYR<=64','Age 65-74'='RIDAGEYR>=65 & RIDAGEYR<=74')
subs<-do.call(rbind,lapply(names(groups),function(n){g<-tryCatch(subset(main,eval(parse(text=groups[[n]]))),error=function(e) NULL);if(is.null(g))return(NULL); r<-rbind(metric(g,'case==1','positive','Sensitivity'),metric(g,'case==0','I(1-positive)','Specificity'));r$group<-n;r$events<-sum(g$variables$case);r}))
write.csv(subs,file.path(out,'subgroup_exploratory.csv'),row.names=FALSE,fileEncoding='UTF-8')
writeLines(c(paste('R:',R.version.string),paste('survey:',packageVersion('survey')),paste('foreign:',packageVersion('foreign')),'AUC: point estimate only; 95% CI unavailable.','Subgroup intervals may be unstable and are exploratory.','Pregnancy status publicly available only for women aged 20-44; 45-74 cannot be classified as not pregnant.'),file.path(out,'analysis_notes.txt'),useBytes=TRUE)
cat('DONE: outputs written to',out,'\n')

