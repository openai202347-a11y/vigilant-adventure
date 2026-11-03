from pathlib import Path
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib import rcParams
rcParams.update({'font.family':'DejaVu Sans','font.size':9,'ps.fonttype':42,'pdf.fonttype':42,'axes.linewidth':0.7})
p=argparse.ArgumentParser();p.add_argument('--paper-dir',default=str(Path(__file__).resolve().parents[1]));p.add_argument('--output-dir',required=True);args=p.parse_args()
r=Path(args.paper_dir);out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
flow=pd.read_csv(r/'04_analysis/results/flow.csv'); metrics=pd.read_csv(r/'04_analysis/results/screening_metrics.csv').set_index('metric')
assert len(flow)==8 and list(flow.remaining)==[11933,6411,5340,4973,3936,3869,3794,3794]
labels=['Interview records','Aged 25–74 years','No self-reported diabetes','Pregnancy criteria met','HbA1c available','BMI available and ≥18 kg/m²','Blood pressure component available','Valid survey weight and design']
fig,ax=plt.subplots(figsize=(6.8,9));ax.set_xlim(0,10.5);ax.set_ylim(-.6,16.3);ax.axis('off')
for i,(label,n,excluded) in enumerate(zip(labels,flow.remaining,flow.excluded)):
 y=15.5-i*2.05;ax.add_patch(FancyBboxPatch((0.2,y-.65),6.5,1.1,boxstyle='round,pad=.08,rounding_size=.10',facecolor='white',edgecolor='black',lw=.8));ax.text(3.45,y+.12,label,ha='center',va='center',fontsize=9);ax.text(3.45,y-.30,f'n = {n:,}',ha='center',va='center',fontsize=9,weight='bold')
 if i<7:ax.annotate('',xy=(3.45,y-1.39),xytext=(3.45,y-.75),arrowprops={'arrowstyle':'->','color':'black','lw':.8})
 if i>0:ax.text(7,y-.09,f'Excluded: {excluded:,}',ha='left',va='center',fontsize=8.5)
fig.savefig(out/'Fig1.eps',bbox_inches='tight');fig.savefig(out/'Fig1_preview.png',dpi=300,bbox_inches='tight');plt.close(fig)
keys=['Sensitivity','Specificity','Positive predictive value','Negative predictive value'];df=metrics.loc[keys];fig,ax=plt.subplots(figsize=(6.9,3.2));y=list(range(4))
ax.errorbar(df.estimate,y,xerr=[df.estimate-df.lower,df.upper-df.estimate],fmt='o',markersize=4,capsize=3,color='black',elinewidth=.9);ax.set_yticks(y,keys);ax.invert_yaxis();ax.set_xlim(0,1.04);ax.set_xlabel('Survey-weighted proportion (95% CI)');ax.grid(axis='x',color='0.85',linewidth=.5);fig.tight_layout()
fig.savefig(out/'Fig2.eps',bbox_inches='tight');fig.savefig(out/'Fig2_preview.png',dpi=300,bbox_inches='tight');plt.close(fig)
print('Created vector EPS and previews from flow.csv and screening_metrics.csv.')
