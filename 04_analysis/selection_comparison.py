from pathlib import Path
import pandas as pd,numpy as np
exec(Path(__file__).with_name('independent_check.py').read_text(encoding='utf-8-sig').split('a=a[keep]')[0])
base=a.between(25,74)&d.DIQ010.eq(2)
assert int(base.sum())==5340 and int(keep.sum())==3794
rows=[]
for label,sel in [('Included',base & keep),('Excluded',base & ~keep)]:
 g=d.loc[sel];n=len(g)
 for name,mask in [('Female',g.RIAGENDR.eq(2)),('Male',g.RIAGENDR.eq(1)),('Age 25-44',g.RIDAGEYR.between(25,44)),('Age 45-64',g.RIDAGEYR.between(45,64)),('Age 65-74',g.RIDAGEYR.between(65,74))]:
  rows.append(dict(group=label,n_group=n,characteristic=name,n=int(mask.sum()),percent=100*mask.sum()/n,mean_years=np.nan))
 rows.append(dict(group=label,n_group=n,characteristic='Age mean',n=n,percent=np.nan,mean_years=g.RIDAGEYR.mean()))
out=pd.DataFrame(rows);out.to_csv(r/'04_analysis/results/selection_comparison_unweighted.csv',index=False);print(out.to_string(index=False))
