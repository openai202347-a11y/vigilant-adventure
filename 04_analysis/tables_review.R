args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2)stop('Usage: Rscript tables_review.R PAPER_DIR R_LIBRARY')
source(file.path(args[1],'04_analysis','run.R'),encoding='UTF-8')
# Display cell counts and survey weighted means/proportions; analytical sample is a survey domain.
z<-main$variables
catcounts<-function(name,grp){
 labels<-sort(unique(grp[!is.na(grp)])); counts<-table(factor(grp,levels=labels))
 tmp<-paste0('cat',gsub('[^A-Za-z0-9]','',name)); main$variables[[tmp]]<<-factor(grp,levels=labels)
 z<-svymean(as.formula(paste0('~',tmp)),main,na.rm=TRUE)
 data.frame(variable=name,category=labels,n=as.integer(counts),weighted_percent=100*as.numeric(z),se_percent=100*as.numeric(SE(z)))
}
tables<-rbind(catcounts('Sex',ifelse(z$RIAGENDR==1,'Male','Female')),catcounts('Age group',ifelse(z$RIDAGEYR<45,'25-44',ifelse(z$RIDAGEYR<65,'45-64','65-74'))),catcounts('NDS class',ifelse(z$positive==1,'28 or higher','below 28')),catcounts('Hypertension',ifelse(z$highbp==1,'Present','Absent')),catcounts('HbA1c',ifelse(z$case==1,'at least 6.5%','below 6.5%')))
write.csv(tables,file.path(out,'participant_characteristics.csv'),row.names=FALSE)
meanage<-svymean(~RIDAGEYR,main);meanbmi<-svymean(~BMXBMI,main)
write.csv(data.frame(variable=c('Age years','BMI kg/m2'),mean=c(coef(meanage),coef(meanbmi)),se=c(SE(meanage),SE(meanbmi))),file.path(out,'participant_means.csv'),row.names=FALSE)
ct<-as.data.frame(table(score_class=ifelse(z$positive==1,'Screen positive','Screen negative'),outcome=ifelse(z$case==1,'HbA1c >=6.5','HbA1c <6.5')))
write.csv(ct,file.path(out,'unweighted_2x2.csv'),row.names=FALSE)
# Approximate effective sample sizes based on design variance are diagnostic only, not complete NCHS reliability classification.
relfun<-function(domain,expr,label){q<-subset(main,eval(parse(text=domain))); f<-as.formula(paste0('~',expr)); p<-svymean(f,q,na.rm=TRUE); v<-as.numeric(SE(p)^2); e<-as.numeric(coef(p)); data.frame(metric=label,n=nrow(q$variables),psus=length(unique(paste(q$variables$SDMVSTRA,q$variables$SDMVPSU))),strata=length(unique(q$variables$SDMVSTRA)),df=degf(q),p=e,se=sqrt(v),n_eff=ifelse(v>0,e*(1-e)/v,NA))}
r<-do.call(rbind,list(relfun('case==1','positive','Sensitivity'),relfun('case==0','I(1-positive)','Specificity'),relfun('positive==1','case','PPV'),relfun('positive==0','I(1-case)','NPV'),relfun('TRUE','case','Prevalence')))
r$n_eff_capped<-pmin(r$n,r$n_eff)
r$n_eff_adjusted<-pmin(r$n,r$n_eff*(qt(.025,r$n-1)/qt(.025,r$df))^2)
r$kg_lower<-qbeta(.025,r$n_eff_adjusted*r$p,r$n_eff_adjusted*(1-r$p)+1)
r$kg_upper<-qbeta(.975,r$n_eff_adjusted*r$p+1,r$n_eff_adjusted*(1-r$p))
r$kg_absolute_width<-r$kg_upper-r$kg_lower
r$kg_relative_width<-r$kg_absolute_width/r$p
r$flag<-ifelse(r$n<30 | r$n_eff_capped<30 | r$df<8 | r$kg_absolute_width>=.30 | (r$kg_absolute_width>.05 & r$kg_relative_width>1.30),'Review/suppress under NCHS criteria','Meets checked numerical criteria; independent review required')
write.csv(r,file.path(out,'reliability_diagnostics.csv'),row.names=FALSE)
cat('Supplementary participant and diagnostic tables created.\n')


