"""Recompute the manuscript summary statistics from results.pkl."""
import pickle, numpy as np, pandas as pd, itertools, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
import validation_engine as ve
RDLogger.DisableLog('rdApp.*')
o=pickle.load(open('./results.pkl','rb'))
M=o['main']; order=list(dict.fromkeys(M.Prototype))
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
fp=lambda s: gen.GetFingerprint(Chem.MolFromSmiles(s))
def intdiv(smis):
    f=[fp(s) for s in smis]; v=[1-DataStructs.TanimotoSimilarity(a,b) for a,b in itertools.combinations(f,2)]
    return float(np.mean(v))
R={}
# --- main headline
print('MAIN records',len(M),'distinct',M.SMILES.nunique(),'direct',(M.EtapaGeracao=='protótipo').sum(),'expl',(M.EtapaGeracao!='protótipo').sum())
print('score %.2f/%.2f tan %.3f/%.3f sa %.2f/%.2f alertfree %d'%(M.ScoreGlobal.mean(),M.ScoreGlobal.std(),M.SimilaridadeTanimoto.mean(),M.SimilaridadeTanimoto.std(),M.SAScore.mean(),M.SAScore.std(),(M.QtdAlertas==0).sum()))
print('cat',M.CategoriaScore.value_counts().to_dict())
# --- Table 3: overlap & diversity
rows=[]
for p in order:
    d=M[M.Prototype==p]; sets=[set(d[d.Seed==s].SMILES) for s in (41242,41243,41244)]
    U=set.union(*sets); I=set.intersection(*sets)
    dr=d[d.EtapaGeracao=='protótipo'].groupby('Seed').size().mean(); ex=d[d.EtapaGeracao!='protótipo'].groupby('Seed').size().reindex([41242,41243,41244],fill_value=0).mean()
    jac=np.mean([len(a&b)/len(a|b) for a,b in itertools.combinations(sets,2)])
    idv=np.mean([intdiv(list(s)) for s in sets])
    scaf=len({MurckoScaffold.MurckoScaffoldSmiles(s) for s in U})
    rows.append(dict(Prototype=p,Direct=dr,Explor=ex,Union=len(U),Shared=len(I),Jaccard=jac,IntDiv=idv,Scaffolds=scaf))
T3=pd.DataFrame(rows); print(T3.round(3).to_string()); R['T3']=T3
print('pooled intdiv mean %.3f'%T3.IntDiv.mean())
# --- funnel
mm=o['main_meta']; mm['acc']=mm.returned/mm.attempts
print(mm.groupby('Prototype').agg(att=('attempts','mean'),pool=('pool','mean'),acc=('acc','mean')).loc[order].round(3))
# --- recompute score to validate Table 1
def comps(r):
    lip=0;veb=0  # all pass hard filters
    dl=0.65*100+0.35*100
    ph=np.mean([ve._range_score(r.MW,180,500,180),ve._range_score(r.LogP,0.0,4.5,2.0),ve._range_score(r.TPSA,20,120,80)])
    sim=min(100,max(0,r.SimilaridadeTanimoto*100))
    sy=ve.EASE_SCORE[str(r.FacilidadeSintese).lower()]
    sa=min(100,max(0,100-(r.SAScore-1)/9*100))
    al=min(100,max(0,100-30*r.QtdAlertas))
    return [dl,ph,sim,sy,sa,al]
C=np.array([comps(r) for r in M.itertuples()])
w0=np.array([.18,.18,.18,.18,.14,.14]); base=C@w0
print('recon max abs diff vs ScoreGlobal: %.3f'%np.abs(base-M.ScoreGlobal.values).max())
print('component means',C.mean(0).round(1),'physchem min',C[:,1].min().round(1),'sd',C.std(0).round(1))
# --- weight sensitivity
rng=np.random.default_rng(20261017)
runs=[g.index.values for _,g in M.groupby(['Prototype','Seed'])]
rho=[];top5=[]
for _ in range(2000):
    w=w0*rng.uniform(0.5,1.5,6); w/=w.sum()
    s=C@w
    for ix in runs:
        a=pd.Series(base[ix]).rank().values; b=pd.Series(s[ix]).rank().values
        rho.append(np.corrcoef(a,b)[0,1])
        top5.append(len(set(np.argsort(-base[ix])[:5])&set(np.argsort(-s[ix])[:5]))/5)
rho=np.array(rho);top5=np.array(top5)
print('Spearman median %.3f IQR %.3f-%.3f p5 %.3f; top5 overlap median %.2f mean %.2f p5 %.2f'%(np.median(rho),*np.percentile(rho,[25,75]),np.percentile(rho,5),np.median(top5),top5.mean(),np.percentile(top5,5)))
R['sens']=dict(rho_med=float(np.median(rho)),rho_q1=float(np.percentile(rho,25)),rho_q3=float(np.percentile(rho,75)),rho_p5=float(np.percentile(rho,5)),t5_mean=float(top5.mean()),t5_med=float(np.median(top5)))
# per-component leave-out (set weight 0, renormalize)
for i,nm in enumerate(['druglike','physchem','similarity','synthesis','SA','alerts']):
    w=w0.copy(); w[i]=0; w/=w.sum(); s=C@w
    r_=[np.corrcoef(pd.Series(base[ix]).rank(),pd.Series(s[ix]).rank())[0,1] for ix in runs]
    print('drop',nm,'median rho %.3f'%np.median(r_))
# --- direct yield & baselines
print('direct yield (target 200):',o['direct_yield'])
U_=o['uniform']
def expl_stats(d):
    e=d[d.EtapaGeracao!='protótipo']
    return dict(n=len(e),distinct=e.SMILES.nunique(),tan=e.SimilaridadeTanimoto.mean(),score=e.ScoreGlobal.mean(),sa=e.SAScore.mean(),af=100*(e.QtdAlertas==0).mean(),t=e.PrioridadeSintese.mean())
a=expl_stats(M[M.Prototype.isin(['Paracetamol','Ibuprofen','Caffeine','Lidocaine','Benzocaine'])]); b=expl_stats(U_[U_.Prototype.isin(['Paracetamol','Ibuprofen','Caffeine','Lidocaine','Benzocaine'])])
print('weighted exploration',{k:round(v,3) for k,v in a.items()}); print('uniform exploration',{k:round(v,3) for k,v in b.items()}); R['expl']=(a,b)
print('uniform total distinct',U_.SMILES.nunique(),'weighted',M.SMILES.nunique(), 'uniform direct',(U_.EtapaGeracao=='protótipo').sum())
# --- seed10
S10=o['seed10']; rows=[]
for p in order:
    d=S10[S10.Prototype==p]; sets=[set(g.SMILES) for _,g in d.groupby('Seed')]
    rows.append(dict(Prototype=p,Union10=len(set.union(*sets)),Shared10=len(set.intersection(*sets)),Jac=np.mean([len(x&y)/len(x|y) for x,y in itertools.combinations(sets,2)])))
T10=pd.DataFrame(rows); print(T10.round(3).to_string()); R['T10']=T10
# --- extended panel
E=o['ext']; em=o['ext_meta']
print('EXT runs',len(em),'all reach N',(em.returned==30).all(),'records',len(E),'distinct',E.SMILES.nunique())
sat=[];
for p,d in E.groupby('Prototype'):
    sets=[set(g.SMILES) for _,g in d.groupby('Seed')]; sat.append((p,len(sets[0]|sets[1]|sets[2])==30, d.ScoreGlobal.mean(),(d.QtdAlertas==0).mean()*100,d.SimilaridadeTanimoto.mean(),(d.EtapaGeracao!='protótipo').mean()*100))
Es=pd.DataFrame(sat,columns=['Prototype','Saturated','Score','AlertFree','Tan','ExplPct']); print(Es.round(2).to_string())
print('ext: saturated',Es.Saturated.sum(),'of',len(Es),'; mean score %.2f sd %.2f; tan %.3f; alertfree %.1f; explshare %.1f'%(E.ScoreGlobal.mean(),E.ScoreGlobal.std(),E.SimilaridadeTanimoto.mean(),100*(E.QtdAlertas==0).mean(),100*(E.EtapaGeracao!='protótipo').mean()))
print('ext categories',E.CategoriaScore.value_counts().to_dict(),'tiers',E.PrioridadeSintese.value_counts().sort_index().to_dict())
R['Es']=Es
pickle.dump(R,open('./stats.pkl','wb'))
