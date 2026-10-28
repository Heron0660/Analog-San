"""Reproducibility script: Analog-San v1.3 evaluation (main benchmark and complementary analyses).
Requires chem_engine and validation_engine importable, RDKit 2026.03.1, Python 3.13."""
import re, pickle, itertools, random, sys, json, hashlib
from pathlib import Path

import pandas as pd
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import chem_engine as ce, validation_engine as ve
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
SEEDS=(41242,41243,41244)
ROUTES=[ce.ROUTE_CLASSIC,ce.ROUTE_NON_CLASSIC]
MAIN={'Paracetamol':'CC(=O)Nc1ccc(O)cc1','Aspirin':'CC(=O)Oc1ccccc1C(=O)O',
'Ibuprofen':'CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O','Caffeine':'Cn1cnc2c1c(=O)n(C)c(=O)n2C',
'Lidocaine':'CCN(CC)CC(=O)Nc1c(C)cccc1C','Benzocaine':'CCOC(=O)c1ccc(N)cc1',
'Propranolol':'CC(C)NCC(O)COc1cccc2ccccc12','Diclofenac':'OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl'}
EXT={'Naproxen':'COc1ccc2cc(ccc2c1)C(C)C(=O)O','Ketoprofen':'CC(C(=O)O)c1cccc(c1)C(=O)c1ccccc1',
'Indomethacin':'COc1ccc2n(C(=O)c3ccc(Cl)cc3)c(C)c(CC(=O)O)c2c1','Salicylic acid':'OC(=O)c1ccccc1O',
'Phenacetin':'CCOc1ccc(NC(C)=O)cc1','Procaine':'CCN(CC)CCOC(=O)c1ccc(N)cc1','Nicotinamide':'NC(=O)c1cccnc1',
'Theophylline':'Cn1c(=O)c2[nH]cnc2n(C)c1=O','Metoprolol':'COCCc1ccc(OCC(O)CNC(C)C)cc1',
'Atenolol':'CC(C)NCC(O)COc1ccc(CC(N)=O)cc1','Warfarin':'CC(=O)CC(c1ccccc1)C1=C(O)c2ccccc2OC1=O',
'Diazepam':'CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc12','Fluoxetine':'CNCCC(Oc1ccc(cc1)C(F)(F)F)c1ccccc1',
'Celecoxib':'Cc1ccc(cc1)-c1cc(nn1-c1ccc(cc1)S(N)(=O)=O)C(F)(F)F','Haloperidol':'OC1(CCN(CCCC(=O)c2ccc(F)cc2)CC1)c1ccc(Cl)cc1',
'Metoclopramide':'CCN(CC)CCNC(=O)c1cc(Cl)c(N)cc1OC','Tramadol':'COc1cccc(c1)C1(O)CCCCC1CN(C)C',
'Phenytoin':'O=C1NC(=O)C(N1)(c1ccccc1)c1ccccc1','Carbamazepine':'NC(=O)N1c2ccccc2C=Cc2ccccc12',
'Ciprofloxacin':'OC(=O)C1=CN(C2CC2)c2cc(N3CCNCC3)c(F)cc2C1=O'}
pat=re.compile(r'Concluído: (\d+)/(\d+) análogos em (\d+) tentativas\s*\|\s*Pool final: (\d+)')
def run(name,smi,seed,n=30,veber=True):
    logs=[]
    res=ce.generate_analogs(smi,n,veber,ROUTES,log_callback=logs.append,random_seed=seed)
    m=[pat.search(l) for l in logs]; m=[x for x in m if x]
    att=int(m[-1].group(3)) if m else None; pool=int(m[-1].group(4)) if m else None
    rows=[dict(r,Prototype=name,Seed=seed) for r in res]
    return rows,dict(Prototype=name,Seed=seed,N=n,returned=len(res),attempts=att,pool=pool)
out={}
rows=[];meta=[]
for nme,smi in MAIN.items():
    for s in SEEDS:
        r,m=run(nme,smi,s); rows+=r; meta.append(m)
out['main']=pd.DataFrame(rows); out['main_meta']=pd.DataFrame(meta)
# extended panel
rows=[];meta=[]
for nme,smi in EXT.items():
    for s in SEEDS:
        r,m=run(nme,smi,s); rows+=r; meta.append(m)
out['ext']=pd.DataFrame(rows); out['ext_meta']=pd.DataFrame(meta)
# 10 seeds for main panel
rows=[]
for nme,smi in MAIN.items():
    for s in range(41242,41252):
        r,m=run(nme,smi,s); rows+=r
out['seed10']=pd.DataFrame(rows)
# direct-stage full yield (target 200)
dy={}
for nme,smi in MAIN.items():
    r,m=run(nme,smi,41242,n=200)
    d=pd.DataFrame(r); dy[nme]=dict(direct=int((d.EtapaGeracao=='protótipo').sum()),total=len(d),attempts=m['attempts'])
out['direct_yield']=dy
# Uniform-sampling baseline: keep direct-stage ordering, tier metadata, and ease
# labels unchanged; only ignore transformation weights during exploration.
class UniformRandom(random.Random):
    def choices(self, population, weights=None, *, cum_weights=None, k=1):
        return super().choices(population, weights=None, cum_weights=None, k=k)

original_random_class = ce.random.Random
ce.random.Random = UniformRandom
rows=[]
try:
    for nme,smi in MAIN.items():
        for s in SEEDS:
            r,m=run(nme,smi,s); rows+=r
finally:
    ce.random.Random = original_random_class
out['uniform']=pd.DataFrame(rows)

study_dir = Path(__file__).resolve().parents[1]
data_dir = study_dir / 'dados'
data_dir.mkdir(exist_ok=True)
out['main'].to_csv(data_dir / 'records_main_benchmark.csv', index=False, encoding='utf-8-sig')
out['ext'].to_csv(data_dir / 'records_extended_panel.csv', index=False, encoding='utf-8-sig')
out['seed10'].to_csv(data_dir / 'records_ten_seeds.csv', index=False, encoding='utf-8-sig')
out['uniform'].to_csv(data_dir / 'records_uniform_sampling_baseline.csv', index=False, encoding='utf-8-sig')
pd.DataFrame.from_dict(out['direct_yield'], orient='index').rename_axis('Prototype').reset_index().to_csv(
    data_dir / 'direct_enumeration_yield.csv', index=False, encoding='utf-8-sig'
)
out['main_meta'].to_csv(data_dir / 'run_metadata_main.csv', index=False, encoding='utf-8-sig')
out['ext_meta'].to_csv(data_dir / 'run_metadata_extended.csv', index=False, encoding='utf-8-sig')

metadata = {
    'software': f'Analog-San {ve.APP_VERSION} (validation module {ve.VALIDATION_VERSION})',
    'python': sys.version.split()[0],
    'rdkit': __import__('rdkit').__version__,
    'pandas': pd.__version__,
    'numpy': np.__version__,
    'main_seeds': list(SEEDS),
    'ten_seeds': list(range(41242, 41252)),
    'target_per_run': 30,
    'routes': ['classical', 'non-classical'],
    'veber_filter': True,
    'baseline': 'UniformRandom overrides only random-choice weights; priority tiers and direct-stage ordering are unchanged.',
    'sha256': {},
}
for path in sorted(data_dir.glob('*.csv')):
    metadata['sha256'][path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
(study_dir / 'metadata_run.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
pickle.dump(out,open('./results.pkl','wb'))
print({k:(v.shape if hasattr(v,'shape') else len(v)) for k,v in out.items()})
