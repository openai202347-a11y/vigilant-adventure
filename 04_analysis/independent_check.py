from pathlib import Path
import pandas as pd, numpy as np
r=Path(__file__).resolve().parents[1]; src=r/'03_data/raw'
def load(n, cols):
 d=pd.read_sas(src/(n+'.xpt'),format='xport',encoding='utf-8'); return d[cols]
d=load('DEMO_L',['SEQN','RIDAGEYR','RIAGENDR','RIDEXPRG','SDMVPSU','SDMVSTRA'])
for n,cols in [('DIQ_L',['DIQ010']),('GHB_L',['LBXGH','WTPH2YR']),('BMX_L',['BMXBMI']),('BPXO_L',['BPXOSY1','BPXOSY2','BPXOSY3','BPXODI1','BPXODI2','BPXODI3']),('BPQ_L',['BPQ020','BPQ150'])]: d=d.merge(load(n,['SEQN']+cols),on='SEQN',how='left',validate='one_to_one')
a=d.RIDAGEYR;b=d.BMXBMI;s=d[['BPXOSY1','BPXOSY2','BPXOSY3']].mean(axis=1);t=d[['BPXODI1','BPXODI2','BPXODI3']].mean(axis=1)
high=(s.ge(140)|t.ge(90)|d.BPQ150.eq(1)); low=s.lt(140)&t.lt(90)&(d.BPQ020.eq(2)|d.BPQ150.eq(2)); bp=high|low
keep=a.between(25,74)&d.DIQ010.eq(2)&(~(d.RIAGENDR.eq(2)&a.between(25,44))|d.RIDEXPRG.eq(2))&d.LBXGH.notna()&b.ge(18)&bp&d.WTPH2YR.gt(0)&d.SDMVPSU.notna()&d.SDMVSTRA.notna()
a=a[keep];b=b[keep];d=d.loc[keep].copy();high=high[keep]
d['score']=np.select([a.lt(35),a.lt(45),a.lt(55),a.lt(65)],[0,8,13,20],22)+np.select([b.lt(20),b.lt(23),b.lt(25),b.lt(27.5),b.lt(30)],[-6,0,9,14,18],20)+8*high.astype(int)
y=d.LBXGH.ge(6.5);p=d.score.ge(28);w=d.WTPH2YR
se=w[y&p].sum()/w[y].sum();sp=w[~y&~p].sum()/w[~y].sum();ppv=w[y&p].sum()/w[p].sum()
rmetrics=pd.read_csv(r/'04_analysis/results/screening_metrics.csv').set_index('metric')
assert len(d)==3794 and y.sum()==93
for k,v in [('Sensitivity',se),('Specificity',sp),('Positive predictive value',ppv)]: assert abs(v-rmetrics.loc[k,'estimate'])<1e-12,(k,v)
print('Independent reconstruction PASS:',len(d),'participants,',y.sum(),'positive; weighted sensitivity, specificity and PPV agree.')
