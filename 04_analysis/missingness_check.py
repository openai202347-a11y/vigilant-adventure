from pathlib import Path
import pandas as pd,numpy as np
src=Path(__file__).with_name('independent_check.py').read_text(encoding='utf-8-sig').split('a=a[keep]')[0];exec(src)
# Sequential exclusions distinguish structural non-response and coded categories. Values are unweighted.
sel=a.between(25,74); rows=[]
def rec(stage,label,mask): rows.append({'population':stage,'category':label,'count':int(mask.fillna(False).sum())})
for val,label in [(1,'Yes'),(2,'No'),(3,'Borderline'),(7,'Refused'),(9,'Do not know')]:rec('Age 25-74','DIQ010: '+label,sel&d.DIQ010.eq(val))
rec('Age 25-74','DIQ010: missing',sel&d.DIQ010.isna())
base=sel&d.DIQ010.eq(2); subgroup=base&d.RIAGENDR.eq(2)&a.between(25,44)
for val,label in [(1,'Pregnant'),(2,'Not pregnant'),(3,'Indeterminate')]:rec('Women 25-44 without diagnosed diabetes','Pregnancy: '+label,subgroup&d.RIDEXPRG.eq(val))
rec('Women 25-44 without diagnosed diabetes','Pregnancy: missing',subgroup&d.RIDEXPRG.isna())
base=base&(~(d.RIAGENDR.eq(2)&a.between(25,44))|d.RIDEXPRG.eq(2))
rec('After pregnancy criteria','HbA1c missing',base&d.LBXGH.isna());base=base&d.LBXGH.notna()
rec('After HbA1c criteria','BMI missing',base&d.BMXBMI.isna());rec('After HbA1c criteria','BMI below 18',base&d.BMXBMI.notna()&d.BMXBMI.lt(18));base=base&d.BMXBMI.ge(18)
rec('After BMI criteria','Hypertension indeterminate',base&~bp)
rec('After BMI criteria','Systolic pressure unavailable',base&s.isna());rec('After BMI criteria','Diastolic pressure unavailable',base&t.isna())
res=pd.DataFrame(rows);out=r/'04_analysis/results/missingness_and_exclusions.csv';res.to_csv(out,index=False);print(res.to_string(index=False))
