#!/usr/bin/env python
"""fast-mic vs SMETANA head-to-head (run in conda env `smetana`).
For a sampled set of probiotic×UHGG pairs, compute SMETANA MIP/MRO on the fast-mic
medium and join with fast-mic C / interaction_type. Compatibility fixes:
 objective->R_bio1; biomass-drain R_EX_cpd11416_c0 -> SINK; medium fmt e0_pool."""
import sys, os, csv, time, random, warnings
warnings.filterwarnings("ignore")
from math import inf
from reframed import load_cbmodel, Environment, set_default_solver
from reframed.core.model import ReactionType
from smetana.legacy import Community
set_default_solver("cplex")   # else reframed auto-picks gurobi>cplex>scip
from smetana.smetana import mip_score, mro_score

FM="/Users/jingyi/fast-mic"
LEVEL=sys.argv[1] if len(sys.argv)>1 else "L5_pectin"
SYS=sys.argv[2] if len(sys.argv)>2 else "lac"
N=int(sys.argv[3]) if len(sys.argv)>3 else 60
MEDCSV={"L5_pectin":"gradient_L5_pectin_gapseq.csv","L6_resistant_starch":"gradient_L6_resistant_starch_gapseq.csv"}[LEVEL]
med=[(r["compounds"], float(r["maxFlux"])) for r in csv.DictReader(open(f"{FM}/media/{MEDCSV}"))]
cpds=[c for c,_ in med]
probdir={"akk":f"{FM}/test/akk/akk_genomes_faa_gapseq_wdm_xml","lac":f"{FM}/test/lac/lac_genomes_faa_gapseq_wdm_xml"}[SYS]
uhggdir=f"{FM}/test/UHGG/final_gapseq_xml"

# sample pairs stratified by interaction_type from fast-mic
rows=[r for r in csv.DictReader(open(f"{FM}/results/{SYS}_vs_uhgg/{LEVEL}.tsv"),delimiter="\t")
      if float(r["growth_a_alone"])>1e-4 and float(r["growth_b_alone"])>1e-4]
random.seed(42)
mut=[r for r in rows if r["interaction_type"]=="mutualism"]
comp=[r for r in rows if r["interaction_type"]=="competition"]
sample=random.sample(mut,min(N//2,len(mut)))+random.sample(comp,min(N//2,len(comp)))
random.shuffle(sample)

_cache={}
def prep(path,mid):
    if mid in _cache: return _cache[mid]
    m=load_cbmodel(path, flavor="fbc2"); m.id=mid
    m.set_objective({"R_bio1":1.0}); m.biomass_reaction="R_bio1"
    if "R_EX_cpd11416_c0" in m.reactions: m.reactions["R_EX_cpd11416_c0"].reaction_type=ReactionType.SINK
    _cache[mid]=m; return m

out=open(f"{FM}/results/smetana_comparison/smetana_{SYS}_{LEVEL}.tsv","w")
w=csv.writer(out,delimiter="\t")
w.writerow(["species_a","species_b","fm_type","fm_C","fm_benefit_a","fm_benefit_b","fm_xfeed","smetana_MIP","smetana_MRO","sec"])
maxup=10.0
done=0
for r in sample:
    a,b=r["species_a"],r["species_b"]
    pa=f"{probdir}/{a}.xml"; pb=f"{uhggdir}/{b}.xml"
    if not(os.path.exists(pa) and os.path.exists(pb)): continue
    t0=time.time()
    try:
        m1=prep(pa,a); m2=prep(pb,b)
        comm=Community("c",[m1,m2], copy_models=True)
        env=Environment()
        for cpd, mf in med:
            env[f"R_EX_M_{cpd}_e0_pool"] = (-mf, inf)   # uptake ≤ fast-mic per-compound maxFlux
        mip=mip_score(comm, environment=env, min_growth=0.01, max_uptake=maxup, verbose=False)
        mro=mro_score(comm, environment=env, min_growth=0.01, max_uptake=maxup, verbose=False)
        MIP=mip[0] if isinstance(mip,tuple) else mip
        MRO=mro[0] if isinstance(mro,tuple) else mro
    except Exception as e:
        MIP=MRO=None
    w.writerow([a,b,r["interaction_type"],r["competition_intensity"],r["benefit_a"],r["benefit_b"],
                r["n_exchanged_metabolites"],MIP,MRO,f"{time.time()-t0:.1f}"]); out.flush()
    done+=1
out.close()
print(f"done {done} pairs -> results/smetana_comparison/smetana_{SYS}_{LEVEL}.tsv")
