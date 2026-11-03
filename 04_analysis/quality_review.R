# Run after the prespecified main analysis, preserving its original outputs.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2) stop('Usage: Rscript quality_review.R PAPER_DIR R_LIBRARY')
source(file.path(args[1], '04_analysis', 'run.R'), encoding='UTF-8')
# Supplemental logit intervals with survey design degrees of freedom.
qmetric <- function(design,domain,expr,label){
 q<-subset(design,eval(parse(text=domain))); deg<-degf(q)
 z<-tryCatch(svyciprop(as.formula(paste0('~',expr)), q, method='logit', df=deg),error=function(e)NULL)
 if(is.null(z)) return(data.frame(metric=label,n=nrow(q$variables),df=deg,estimate=NA,lower=NA,upper=NA))
 ci<-confint(z);data.frame(metric=label,n=nrow(q$variables),df=deg,estimate=unname(coef(z)),lower=ci[1],upper=ci[2])
}
qc<-do.call(rbind,list(qmetric(main,'case==1','positive','Sensitivity'),qmetric(main,'case==0','I(1-positive)','Specificity'),qmetric(main,'positive==1','case','Positive predictive value'),qmetric(main,'positive==0','I(1-case)','Negative predictive value'),qmetric(main,'TRUE','positive','Screen-positive proportion'),qmetric(main,'TRUE','case','Elevated HbA1c prevalence')))
write.csv(qc,file.path(out,'metrics_logit_sensitivity.csv'),row.names=FALSE)
# Supplemental Rao-Wu style bootstrap replicates at PSU-within-strata level.
# Approximation: replicate method is not NHANES official replicate weights; CDC recommends Taylor linearization.
set.seed(20261031)
reps<-as.svrepdesign(des,type='bootstrap',replicates=2000,mse=TRUE)
subrep<-subset(reps,eligible)
xx<-subrep$variables
calc_auc<-function(w){
 v<-xx$score; c<-xx$case; vals<-sort(unique(v))
 wp<-vapply(vals,function(k)sum(w[c==1 & v==k]),numeric(1)); wn<-vapply(vals,function(k)sum(w[c==0 & v==k]),numeric(1))
 if(sum(wp)==0 || sum(wn)==0)return(NA_real_)
 sum(wp*(cumsum(wn)-0.5*wn))/(sum(wp)*sum(wn))
}
basew<-weights(subrep,type='sampling'); br<-weights(subrep,type='analysis')
base_auc<-calc_auc(basew); estimates<-apply(br,2,calc_auc)
if(any(!is.finite(estimates)))stop('At least one replicate lacks both outcome groups')
stopifnot(abs(base_auc-weighted_auc)<1e-10)
# mse=TRUE: sum squared deviation from full estimate scaled by survey replicate scale/rscales.
variance<-reps$scale*sum(reps$rscales*(estimates-base_auc)^2)
se<-sqrt(variance); df<-degf(subrep); crit<-qt(.975,df=df)
auc_out<-data.frame(estimate=base_auc,se=se,lower=max(0,base_auc-crit*se),upper=min(1,base_auc+crit*se),df=df,replicates=length(estimates),seed=20261031,method='Stratified PSU rescaled bootstrap, normal-on-AUC CI using t critical value; supplemental approximation')
write.csv(auc_out,file.path(out,'auc_bootstrap_supplement.csv'),row.names=FALSE)
cat('Supplementary quality review complete; df=',df,' AUC=',base_auc,' SE=',se,'\n')
