args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2)stop('Usage: Rscript auc_linearization_survey_check.R PAPER_DIR R_LIBRARY')
source(file.path(args[1],'04_analysis','run.R'),encoding='UTF-8')
x<-main$variables;w<-as.numeric(weights(main));y<-x$case==1;v<-x$score;lev<-sort(unique(v));W1<-sum(w[y]);W0<-sum(w[!y]);
w1<-vapply(lev,function(t)sum(w[y & v==t]),numeric(1));w0<-vapply(lev,function(t)sum(w[!y & v==t]),numeric(1));
f0<-(cumsum(w0)-0.5*w0)/W0;g1<-(rev(cumsum(rev(w1)))-0.5*w1)/W1;auc<-sum(w1*f0)/W1;idx<-match(v,lev)
influence<-ifelse(y,(f0[idx]-auc)/W1,(g1[idx]-auc)/W0)
main$variables$auc_influence<-influence
z<-svytotal(~auc_influence,main);se<-as.numeric(SE(z)); py<-read.csv(file.path(out,'auc_linearization_exploratory.csv'))
stopifnot(abs(auc-py$estimate)<1e-10,abs(se-py$se)<1e-10)
result<-data.frame(auc=auc,se_R_survey=se,se_python=py$se,df=degf(main),note='Post hoc computational cross-check; no independent human review')
write.csv(result,file.path(out,'auc_linearization_R_crosscheck.csv'),row.names=FALSE)
print(result)
