from pathlib import Path
import numpy as np, pandas as pd
exec(Path(__file__).with_name('check_survey_variance.py').read_text(encoding='utf-8-sig'))
y=d.LBXGH.ge(6.5).to_numpy(); score=d.score.to_numpy();w=d.WTPH2YR.to_numpy(); levels=np.sort(np.unique(score));
w1=np.array([w[y & (score==x)].sum() for x in levels]);w0=np.array([w[(~y) & (score==x)].sum() for x in levels]); W1=w1.sum(); W0=w0.sum(); assert W1>0 and W0>0
f0=(np.cumsum(w0)-0.5*w0)/W0
g1=(np.cumsum(w1[::-1])[::-1]-0.5*w1)/W1
auc=np.sum(w1*f0)/W1
loc=np.searchsorted(levels,score)
influence=np.where(y,w*(f0[loc]-auc)/W1,w*(g1[loc]-auc)/W0)
# Numerical derivative checks the contribution formula for one case and one noncase.
def empirical_auc(ww):
    aw=np.array([ww[y & (score==x)].sum() for x in levels]);bw=np.array([ww[(~y) & (score==x)].sum() for x in levels]);return np.sum(aw*(np.cumsum(bw)-0.5*bw))/(aw.sum()*bw.sum())
for i in [np.flatnonzero(y)[0],np.flatnonzero(~y)[0]]:
    eps=1e-5; pert=w.copy();pert[i]*=(1+eps);derivative=(empirical_auc(pert)-auc)/eps
    assert abs(derivative-influence[i])<1e-8,(i,derivative,influence[i])
tmp=d[['SDMVSTRA','SDMVPSU']].copy();tmp['linear']=influence
sums=tmp.groupby(['SDMVSTRA','SDMVPSU'],as_index=False).linear.sum();clusters=psus.merge(sums,on=['SDMVSTRA','SDMVPSU'],how='left').fillna({'linear':0.0})
variance=0.0
for _,g in clusters.groupby('SDMVSTRA'):
 vals=g.linear.to_numpy();n=len(vals);variance+=n/(n-1)*np.sum((vals-np.mean(vals))**2)
se=np.sqrt(variance)
# No confidence interval is presented: df choice and applicability need specialist review.
reference_auc=pd.read_csv(r/'04_analysis/results/auc_bootstrap_supplement.csv').iloc[0]
assert abs(auc-reference_auc['estimate'])<1e-10
out=pd.DataFrame([dict(method='First-order influence linearization by public strata and PSU; exploratory',estimate=auc,se=se,bootstrap_se=reference_auc['se'])]);out.to_csv(r/'04_analysis/results/auc_linearization_exploratory.csv',index=False)
print(out.to_string(index=False));print('AUC cross-check only; this is exploratory and does not validate bootstrap or replace specialist review.')
