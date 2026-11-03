from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
root=Path('work/sci-writing-harness/papers/2026-10-31-diabetes-screening'); df=pd.read_csv(root/'04_analysis/results/flow.csv');
labels=['All interview records','Aged 25–74 years','No self-reported diabetes','Pregnancy criteria applied','HbA1c available','BMI available and ≥18 kg/m²','Blood pressure component available','Valid survey weights and design']
fig,ax=plt.subplots(figsize=(8,10));ax.set_xlim(0,10);ax.set_ylim(-.3,17.5);ax.axis('off')
for i,(label,n,excluded) in enumerate(zip(labels,df.remaining,df.excluded)):
 y=16-i*2.05; box=FancyBboxPatch((0.5,y-.65),6.5,1.1,boxstyle='round,pad=.10,rounding_size=.12',facecolor='#f5f8fa',edgecolor='#23576b',lw=1.1);ax.add_patch(box); ax.text(3.75,y+.05,label,ha='center',va='center',fontsize=10);ax.text(3.75,y-.32,f'n = {n:,}',ha='center',va='center',fontsize=9,weight='bold');
 if i<7: ax.annotate('',xy=(3.75,y-1.52),xytext=(3.75,y-.76),arrowprops={'arrowstyle':'->','color':'#23576b','lw':1})
 if i>0: ax.text(7.4,y-.08,f'Excluded: {excluded:,}',ha='left',va='center',fontsize=9,color='#4b5563')
ax.set_title('Study population selection',loc='left',fontsize=15,weight='bold',pad=10)
fig.text(.09,.025,'Counts are unweighted; exclusions are sequential.',fontsize=9,color='#4b5563');fig.savefig('outputs/研究样本筛选流程图.png',dpi=300,bbox_inches='tight');plt.close(fig)
