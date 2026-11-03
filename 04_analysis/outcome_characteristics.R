args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2) stop('Usage: Rscript outcome_characteristics.R PAPER_DIR R_LIBRARY')
source(file.path(args[1],'04_analysis','run.R'),encoding='UTF-8')
# Exploratory, post hoc descriptive table. No hypotheses, p-values, or clinical claims.
main<-update(main, sex_female=as.numeric(RIAGENDR==2), age_25_44=as.numeric(RIDAGEYR>=25 & RIDAGEYR<45), age_45_64=as.numeric(RIDAGEYR>=45 & RIDAGEYR<65), age_65_74=as.numeric(RIDAGEYR>=65 & RIDAGEYR<=74), highbp_ind=as.numeric(highbp==1), screen_ind=as.numeric(positive==1))
rows<-list();k<-0
for (v in c(0,1)) {
 q<-subset(main,case==v);n<-nrow(q$variables);label<-if(v==1) 'HbA1c >=6.5%' else 'HbA1c <6.5%'
 for(nm in c('sex_female','age_25_44','age_45_64','age_65_74','highbp_ind','screen_ind')) {
  m<-svymean(as.formula(paste0('~',nm)),q); cell<-sum(q$variables[[nm]]==1);k<-k+1
  rows[[k]]<-data.frame(outcome=label,n_group=n,characteristic=nm,n_cell=cell,weighted_percent=100*as.numeric(coef(m)),se_percent=100*as.numeric(SE(m)),mean=NA,mean_se=NA,flag=ifelse(cell<30 || n-cell<30,'Small unweighted cell; unstable','NCHS reliability not fully assessed'))
 }
 for(nm in c('RIDAGEYR','BMXBMI')) {
  m<-svymean(as.formula(paste0('~',nm)),q);k<-k+1
  rows[[k]]<-data.frame(outcome=label,n_group=n,characteristic=nm,n_cell=n,weighted_percent=NA,se_percent=NA,mean=as.numeric(coef(m)),mean_se=as.numeric(SE(m)),flag='Descriptive; reliability not fully assessed')
 }
}
z<-do.call(rbind,rows);stopifnot(all(z$n_group[z$outcome=='HbA1c >=6.5%']==93),all(z$n_group[z$outcome=='HbA1c <6.5%']==3701))
write.csv(z,file.path(out,'characteristics_by_outcome_exploratory.csv'),row.names=FALSE,na='')
print(z)
