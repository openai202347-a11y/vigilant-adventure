# Independent implementation of stratified-cluster Taylor variance for the pre-existing survey metrics.
from pathlib import Path
exec(Path(__file__).with_name('independent_check.py').read_text(encoding='utf-8-sig'))
full=load('DEMO_L',['SEQN','SDMVSTRA','SDMVPSU']).merge(load('GHB_L',['SEQN','WTPH2YR']),on='SEQN',how='left',validate='one_to_one')
psus=full.loc[full.WTPH2YR.gt(0)&full.SDMVSTRA.notna()&full.SDMVPSU.notna(),['SDMVSTRA','SDMVPSU']].drop_duplicates()
assert psus.groupby('SDMVSTRA').size().min()>1
cases=d.LBXGH.ge(6.5); positive=d.score.ge(28)
checks=[('Sensitivity',cases,positive.astype(float)),('Specificity',~cases,(~positive).astype(float)),('PPV',positive,cases.astype(float)),('NPV',~positive,(~cases).astype(float)),('Prevalence',pd.Series(True,index=d.index),cases.astype(float))]
reference=pd.read_csv(r/'04_analysis/results/reliability_diagnostics.csv').set_index('metric')
for name,domain,indicator in checks:
 sub=d.loc[domain].copy(); obs=indicator.loc[domain]; w=sub.WTPH2YR; estimate=np.sum(w*obs)/np.sum(w)
 sub['linear']=w*(obs-estimate)/np.sum(w)
 sums=sub.groupby(['SDMVSTRA','SDMVPSU'],as_index=False).linear.sum()
 clusters=psus.merge(sums,on=['SDMVSTRA','SDMVPSU'],how='left').fillna({'linear':0.0})
 V=0.0
 for _,g in clusters.groupby('SDMVSTRA'):
  vals=g.linear.to_numpy();n=len(vals);V+=n/(n-1)*np.sum((vals-np.mean(vals))**2)
 standard_error=np.sqrt(V); expected=reference.loc[name,'se']
 print(f'{name}: independent SE={standard_error:.12f}; R survey SE={expected:.12f}; difference={standard_error-expected:.2g}')
 assert abs(standard_error-expected)<1e-10
print('PASS: Taylor-linearised design standard errors agree to 1e-10.')
