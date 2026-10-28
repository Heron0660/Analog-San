"""Generate the five Portuguese manuscript figures from results.pkl."""
import pickle, numpy as np, pandas as pd, matplotlib, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.ticker import FuncFormatter
from rdkit import Chem
from rdkit.Chem import Draw
import chem_engine as ce
OUT='./figuras/'
o=pickle.load(open('./results.pkl','rb')); M=o['main']
order=['Paracetamol','Aspirin','Ibuprofen','Caffeine','Lidocaine','Benzocaine','Propranolol','Diclofenac']
pt=dict(zip(order,['Paracetamol','Aspirina','Ibuprofeno','Cafeína','Lidocaína','Benzocaína','Propranolol','Diclofenaco']))
col=dict(zip(order,['#7FA6CC','#E6C48C','#93BFA9','#D9A9B3','#DDA68B','#B9AED6','#A3ADB8','#C4C892']))
def dark(c,f=0.62):
    c=c.lstrip('#'); r,g,b=[int(c[i:i+2],16) for i in (0,2,4)]
    return '#%02x%02x%02x'%(int(r*f),int(g*f),int(b*f))
mk=dict(zip(order,['o','s','^','D','v','P','X','h']))
INK='#3b4450'
SM={'Paracetamol':'CC(=O)Nc1ccc(O)cc1','Aspirin':'CC(=O)Oc1ccccc1C(=O)O','Ibuprofen':'CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O','Caffeine':'Cn1cnc2c1c(=O)n(C)c(=O)n2C','Lidocaine':'CCN(CC)CC(=O)Nc1c(C)cccc1C','Benzocaine':'CCOC(=O)c1ccc(N)cc1','Propranolol':'CC(C)NCC(O)COc1cccc2ccccc12','Diclofenac':'OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
comma=FuncFormatter(lambda v,p: ('%g'%round(v,6)).replace('.',',').replace('-','−') if abs(v-round(v))>1e-9 else ('%d'%round(v)).replace('-','−'))
def save(fig,n):
    fig.savefig(OUT+n+'.png',dpi=300,bbox_inches='tight'); fig.savefig(OUT+n+'.svg',bbox_inches='tight'); plt.close(fig)
# Fig 1 fluxo
fig,ax=plt.subplots(figsize=(7.4,2.3)); ax.set_xlim(-1,101); ax.set_ylim(0,32); ax.axis('off')
steps=[('Entrada','protótipo\nem SMILES'),('Enumeração','14 transformações\npredefinidas'),('Sanitização','valência (RDKit)\ne SMILES\ncanônico'),('Filtragem','Lipinski;\ncritérios de\nVeber opcionais'),('Priorização','nível da regra\n+ escore\nheurístico'),('Exportação','CSV / SMI / MOL2\n+ figuras do\nrelatório')]
w=14.8; gap=2.1
for i,(t,s_) in enumerate(steps):
    x=i*(w+gap)+0.4
    ax.add_patch(FancyBboxPatch((x,8),w,20,boxstyle='round,pad=0.2,rounding_size=1.2',fc='#EEF3F8',ec='#7C9FC2',lw=1.2))
    ax.text(x+w/2,24,t,ha='center',va='center',fontsize=8,fontweight='bold',color=INK)
    ax.text(x+w/2,15.5,s_,ha='center',va='center',fontsize=6.2,color=INK,linespacing=1.35)
    if i<5: ax.add_patch(FancyArrowPatch((x+w+0.3,18),(x+w+gap-0.3,18),arrowstyle='-|>',mutation_scale=8,color='#7C9FC2',lw=1.2))
ax.text(50,2.5,'Etapa de exploração: seleção aleatória semeada do parente e amostragem de transformações ponderada por prioridade',ha='center',fontsize=7,color='#5b6573')
save(fig,'figura_1_fluxo_de_trabalho')
# Fig 2 escores
fig,ax=plt.subplots(figsize=(7.2,3.8))
data=[M[M.Prototype==p].ScoreGlobal.values for p in order]
bp=ax.boxplot(data,patch_artist=True,whis=1.5,widths=0.55,flierprops=dict(marker='o',markersize=3,markerfacecolor='none',markeredgecolor='#7a828c',markeredgewidth=0.8),medianprops=dict(color=INK,lw=1.6),whiskerprops=dict(color='#7a828c',lw=1),capprops=dict(color='#7a828c',lw=1))
for b_,p in zip(bp['boxes'],order): b_.set_facecolor(col[p]); b_.set_edgecolor(dark(col[p])); b_.set_linewidth(1)
for i,p in enumerate(order,1):
    rm=M[M.Prototype==p].groupby('Seed').ScoreGlobal.mean().values
    ax.scatter(np.full(len(rm),i)+np.linspace(-0.12,0.12,len(rm)),rm,marker='D',s=16,color=INK,zorder=5)
ax.set_xticks(range(1,9)); ax.set_xticklabels([pt[p] for p in order],rotation=25,ha='right'); ax.set_ylabel('ScoreGlobal heurístico (0–100)'); ax.set_xlabel('Protótipo'); ax.set_ylim(60,100)
ax.yaxis.set_major_formatter(comma); ax.yaxis.grid(True,color='#e6e9ed',lw=0.7); ax.set_axisbelow(True)
for sp in ('left','bottom'): ax.spines[sp].set_color('#8a929c')
save(fig,'figura_2_distribuicao_escores')
# Fig 3 MM x MolLogP
fig,ax=plt.subplots(figsize=(7.2,4.6))
for p in order:
    d=M[M.Prototype==p].drop_duplicates('SMILES')
    ax.scatter(d.MW,d.LogP,marker=mk[p],s=30,color=col[p],edgecolor=dark(col[p]),linewidth=0.6,label=pt[p])
for p in order:
    pr=ce.calculate_properties(Chem.MolFromSmiles(SM[p]))
    ax.scatter(pr['MW'],pr['LogP'],marker='*',s=190,color=col[p],edgecolor=INK,linewidth=0.9,zorder=6)
ax.scatter([],[],marker='*',s=130,color='white',edgecolor=INK,label='Composto-mãe (protótipo)')
ax.set_xlabel('Massa molecular (Da)'); ax.set_ylabel('MolLogP de Crippen (RDKit)'); ax.grid(True,color='#e6e9ed',lw=0.7); ax.set_axisbelow(True)
ax.xaxis.set_major_formatter(comma); ax.yaxis.set_major_formatter(comma)
for sp in ('left','bottom'): ax.spines[sp].set_color('#8a929c')
ax.legend(loc='center left',bbox_to_anchor=(1.01,0.5),frameon=False,fontsize=8)
save(fig,'figura_3_mm_mollogp')
# Fig 4 niveis
tc=['#DCE6F0','#BED2E4','#9DB9D4','#7C9FC2','#5C82AB','#3F6590']
fig,ax=plt.subplots(figsize=(7.2,3.9)); bot=np.zeros(8)
ct=pd.crosstab(M.Prototype,M.PrioridadeSintese).reindex(order).reindex(columns=range(1,7),fill_value=0); ct=ct.div(ct.sum(axis=1),axis=0)
for t in range(1,7):
    ax.bar(range(8),ct[t].values,bottom=bot,color=tc[t-1],width=0.7,label=f'Nível {t}',edgecolor='white',linewidth=0.8); bot+=ct[t].values
ax.set_xticks(range(8)); ax.set_xticklabels([pt[p] for p in order],rotation=25,ha='right'); ax.set_ylim(0,1); ax.set_ylabel('Fração de análogos retornados'); ax.set_xlabel('Protótipo')
ax.yaxis.set_major_formatter(FuncFormatter(lambda v,p:('%.1f'%v).replace('.',',')))
for sp in ('left','bottom'): ax.spines[sp].set_color('#8a929c')
ax.legend(title='Nível de modificação interno',ncol=6,loc='lower center',bbox_to_anchor=(0.5,1.01),frameon=False,fontsize=8,title_fontsize=9)
save(fig,'figura_4_niveis_prioridade')
# Fig 5 edicoes
fig,axs=plt.subplots(1,3,figsize=(7.2,2.5))
for ax,smi,t in zip(axs,['CC(=O)Nc1ccc(O)cc1','CC(=O)N(C)c1ccc(O)cc1','COc1ccc(NC(C)=O)cc1'],['Paracetamol (protótipo)','Análogo N-metil','Análogo O-metil']):
    im=Draw.MolToImage(Chem.MolFromSmiles(smi),size=(600,420)); ax.imshow(im); ax.axis('off'); ax.set_title(t,fontsize=9,color=INK,pad=2)
fig.subplots_adjust(wspace=0.02)
save(fig,'figura_5_edicoes_moleculares')
print('ok')
