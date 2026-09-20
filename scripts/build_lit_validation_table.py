#!/usr/bin/env python
"""Literature cross-check of fast-mic's predicted cross-feeding currencies.

Joins the predicted export prevalence (results/fig5/fig5_prevalence_{sys}.tsv, produced by
extract_fig5_crossfeed.py) with hand-curated experimental evidence from the literature, so the
numbers stay reproducible while the evidence column stays auditable.

Verdicts apply to the CURRENCY, i.e. is this metabolite documented as a cross-fed compound in
gut (or gut-like) communities:
  supported     — documented by culture/coculture or gnotobiotic experiment
  partial       — only one half of the exchange (production or consumption) demonstrated
  untested      — no direct experimental report found for this compound as a cross-fed currency
  contradicted  — published experiments point the opposite way to the prediction

Output: results/litvalidation/lit_validation.tsv
Usage:  python3 scripts/build_lit_validation_table.py
"""
import csv, os, collections

FM = "."
OUT = f"{FM}/results/litvalidation/lit_validation.tsv"

# cpd -> (verdict, evidence sentence, organisms/system, method, reference)
EVIDENCE = {
 "cpd00029": ("supported",
   "mucin fermentation by A. muciniphila liberates acetate that becomes directly available to neighbouring "
   "microorganisms, and acetate-based cross-feeding sustains butyrate producers in defined cocultures",
   "A. muciniphila x Anaerostipes caccae / Eubacterium hallii / Faecalibacterium prausnitzii; lactobacilli and "
   "bifidobacteria x colon bacteria",
   "anaerobic coculture, metabolite quantification",
   "Belzer et al. 2017 mBio 8:e00770-17 (10.1128/mbio.00770-17); Moens et al. 2017 Int J Food Microbiol 241:225-236 (10.1016/j.ijfoodmicro.2016.10.019)"),
 "cpd01861": ("supported",
   "A. muciniphila grown on mucin supports syntrophic growth of non-mucolytic butyrate producers; "
   "cobalamin supplied by E. hallii drives A. muciniphila propionate formation via the 1,2-propanediol route",
   "A. muciniphila x Anaerostipes caccae / Eubacterium hallii / Faecalibacterium prausnitzii",
   "anaerobic coculture, metabolite quantification",
   "Belzer et al. 2017 mBio 8:e00770-17 (10.1128/mbio.00770-17); Engels et al. 2016 Front Microbiol 7:713 (10.3389/fmicb.2016.00713)"),
 "cpd00453": ("supported", "see (R)-1,2-propanediol", "as above", "as above",
   "Belzer et al. 2017 mBio (10.1128/mbio.00770-17)"),
 "cpd00116": ("supported",
   "pectin-fermenting Lachnospira multiparus releases methanol; the methylotroph Eubacterium limosum consumes it, "
   "nearly doubling coculture yield over the pectin degrader alone",
   "Lachnospira multiparus x Eubacterium limosum on 0.2% pectin",
   "defined coculture, end-product analysis",
   "Rode, Genthner & Bryant 1981 Appl Environ Microbiol 42:20-22 (10.1128/aem.42.1.20-22.1981)"),
 "cpd00047": ("supported",
   "formate released by the colonic starch degrader Ruminococcus bromii is consumed by acetogens in defined "
   "co-culture, with the exchange confirmed transcriptomically; formate is also a major product of pectin "
   "fermentation in gut isolates",
   "Ruminococcus bromii x Blautia hydrogenotrophica; Lachnospira multiparus on pectin",
   "defined co-culture with transcriptomics; defined co-culture, end-product analysis",
   "Laverde Gomez et al. 2019 Environ Microbiol 21:259-271 (10.1111/1462-2920.14454); Rode et al. 1981 AEM 42:20-22 (10.1128/aem.42.1.20-22.1981)"),
 "cpd00159": ("supported",
   "lactate produced by lactobacilli degrading inulin-type fructans is converted to butyrate by Anaerostipes caccae "
   "in biculture; fecal isolates of E. hallii and A. caccae consume lactate and produce butyrate",
   "Lactobacillus acidophilus x Anaerostipes caccae; human fecal isolates",
   "biculture on oligofructose; isolate characterization",
   "Moens, Verce & De Vuyst 2017 Int J Food Microbiol 241:225-236 (10.1016/j.ijfoodmicro.2016.10.019); Duncan, Louis & Flint 2004 AEM 70:5810-5817 (10.1128/aem.70.10.5810-5817.2004)"),
 "cpd00221": ("supported", "see L-lactate; both enantiomers are used by colonic lactate utilizers",
   "as above", "as above",
   "Duncan et al. 2004 AEM (10.1128/aem.70.10.5810-5817.2004)"),
 "cpd00036": ("supported",
   "succinate released by gut fermenters is decarboxylated to propionate by succinate-utilizing taxa "
   "(Phascolarctobacterium, Dialister, Bacteroides)",
   "human colonic microbiota", "review of propionate/butyrate routes with culture evidence",
   "Louis & Flint 2017 Environ Microbiol 19:29-41 (10.1111/1462-2920.13589)"),
 "cpd00105": ("supported",
   "Bacteroides thetaiotaomicron carries dedicated ribose-utilization systems that scavenge ribose, including "
   "ribose released from nucleosides, and these confer a diet-specific colonization advantage in vivo",
   "B. thetaiotaomicron in gnotobiotic mice",
   "molecular genetics + competitive colonization",
   "Glowacki et al. 2020 Cell Host Microbe 27:79-92 (10.1016/j.chom.2019.11.009)"),
 "cpd00082": ("supported",
   "lactobacilli degrading inulin-type fructans release free fructose into the medium, which becomes available "
   "to non-degrading community members",
   "Lactobacillus acidophilus / L. paracasei on oligofructose and inulin",
   "monoculture and biculture, HPLC of released sugars",
   "Moens et al. 2017 Int J Food Microbiol (10.1016/j.ijfoodmicro.2016.10.019)"),
 "cpd00644": ("supported",
   "B-vitamin auxotrophy is widespread among butyrate producers and is relieved by prototrophic partners in "
   "synthetic cocultures, establishing B vitamins as cross-fed currencies",
   "15 human butyrate-producing strains, synthetic cocultures",
   "genome analysis + defined-medium growth and coculture",
   "Soto-Martin et al. 2020 mBio 11:e00886-20 (10.1128/mbio.00886-20)"),
 "cpd00363": ("supported",
   "both halves of the exchange are documented: ethanol is a major end product of pectin fermentation by gut "
   "isolates, and bacteria of the normal human colonic flora oxidize ethanol, producing acetate",
   "Lachnospira multiparus on pectin; aerobic isolates of the normal human large-intestinal flora",
   "defined co-culture, end-product analysis; isolate assays of alcohol and aldehyde dehydrogenase activity",
   "Rode et al. 1981 AEM 42:20-22 (10.1128/aem.42.1.20-22.1981); Nosova et al. 1996 Alcohol Alcohol 31:555-564 (10.1093/oxfordjournals.alcalc.a008191)"),
 "cpd00118": ("partial",
   "putrescine is produced in the gut and its formation rises on pectin, but the demonstrated producers are "
   "Bacteroides and Fusobacterium rather than the probiotic genera modelled here",
   "B. thetaiotaomicron + F. varium, gnotobiotic rats +/- 10% pectin",
   "gnotobiotic in vivo, cecal polyamine quantification",
   "Noack et al. 2000 J Nutr 130:1225-1231 (10.1093/jn/130.5.1225)"),
 "cpd00033": ("partial",
   "amino-acid auxotrophy and relief by prototrophic partners is documented for gut anaerobes, but exchange of "
   "this specific amino acid was not individually tested",
   "human butyrate producers, synthetic cocultures", "defined-medium growth and coculture",
   "Soto-Martin et al. 2020 mBio (10.1128/mbio.00886-20)"),
 "cpd00053": ("partial", "see glycine", "as above", "as above",
   "Soto-Martin et al. 2020 mBio (10.1128/mbio.00886-20)"),
 "cpd00324": ("untested",
   "colonic bacteria release methanethiol in quantity and the caecal mucosa detoxifies it, so the compound is "
   "certainly produced in the gut; we found no experiment showing another gut bacterium taking it up, which "
   "leaves the exchange itself untested",
   "human and rat caecal mucosa; colonic bacterial flora",
   "gas measurement and mucosal detoxification assays",
   "Levitt et al. 1999 J Clin Invest 104:1107-1114 (10.1172/jci7712)"),
 "cpd00130": ("partial",
   "the consumer side is demonstrated in vivo: Salmonella Typhimurium imports microbiota-derived malate through "
   "DcuABC to drive fumarate respiration during gut-lumen colonization, and lactic acid bacteria carry dedicated "
   "malate transporters; release of malate by a gut partner has not been shown directly",
   "S. Typhimurium in colonized mice; Lactococcus lactis membrane vesicles",
   "in vivo colonization with transporter mutants; transporter biochemistry",
   "Nguyen et al. 2020 Cell Host Microbe 27:922-936 (10.1016/j.chom.2020.04.013); Bandell et al. 1997 J Biol Chem 272:18140-18146 (10.1074/jbc.272.29.18140)"),
 "cpd00020": ("untested",
   "pyruvate is a central intermediate and is set aside in this study as an ambiguous end product; "
   "no direct experimental report of pyruvate cross-feeding in gut communities was located",
   "-", "-", "no direct evidence located (Europe PMC, Sep 2026)"),
 "cpd00246": ("contradicted",
   "published experiments show lactobacilli CONSUMING purine nucleosides rather than exporting them: "
   "Lactiplantibacillus plantarum degrades nucleosides in the gut to lower urate, and Latilactobacillus sakei "
   "uses the pentose moiety of inosine as an energy source. Microbiome-derived inosine is documented, but is "
   "attributed to Bifidobacterium pseudolongum",
   "L. plantarum (mice); L. sakei (in vitro); B. pseudolongum (mice)",
   "in vivo supplementation; growth on nucleosides; gnotobiotic immunotherapy model",
   "Li et al. 2023 Microbiome 11:153 (10.1186/s40168-023-01605-y); Rimaux et al. 2011 AEM 77:6539-6550 (10.1128/aem.00498-11); Mager et al. 2020 Science 369:1481-1489 (10.1126/science.abc3421)"),
}

# Whether the experiment used the probiotic genus modelled here as the DONOR. Kept separate from
# the verdict: a currency can be well documented while the published donor is a different taxon.
DONOR_MATCH = {
 ("Akkermansia", "cpd00029"): "yes",
 ("Akkermansia", "cpd01861"): "yes", ("Akkermansia", "cpd00453"): "yes",
 ("Lactobacillus", "cpd00159"): "yes", ("Lactobacillus", "cpd00221"): "yes",
 ("Lactobacillus", "cpd00082"): "yes",
 ("Akkermansia", "cpd00105"): "consumer side only", ("Lactobacillus", "cpd00105"): "consumer side only",
 ("Akkermansia", "cpd00130"): "consumer side only", ("Lactobacillus", "cpd00130"): "consumer side only",
}

NAME_FIX = {"cpd00324": "Methanethiol", "cpd00644": "Pantothenate (vitamin B5)"}

def prevalence(sys_key):
    agg = collections.defaultdict(list); meta = {}
    with open(f"{FM}/results/fig5/fig5_prevalence_{sys_key}.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            agg[r["cpd"]].append(float(r["prevalence"])); meta[r["cpd"]] = (r["name"], r["class"])
    return agg, meta

def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    rows = []
    for sys_key, label in (("akk", "Akkermansia"), ("lac", "Lactobacillus")):
        agg, meta = prevalence(sys_key)
        for cpd, ev in EVIDENCE.items():
            if cpd not in agg:
                continue
            prevs = agg[cpd]; name, cls = meta[cpd]
            name = NAME_FIX.get(cpd, name)
            if sum(prevs) / len(prevs) < 0.05:
                continue
            rows.append([label, cpd, name, cls,
                         f"{sum(prevs)/len(prevs):.1f}", f"{max(prevs):.1f}", ev[0],
                         DONOR_MATCH.get((label, cpd), "no"), *ev[1:]])
    rows.sort(key=lambda r: (r[0], -float(r[4])))
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["system", "cpd", "metabolite", "class", "mean_export_prevalence_pct",
                    "max_export_prevalence_pct", "verdict", "predicted_donor_taxon_tested",
                    "experimental_evidence", "organisms_tested", "method", "reference"])
        w.writerows(rows)
    n = collections.Counter(r[6] for r in rows)
    print(f"wrote {len(rows)} rows -> {OUT}")
    print("  " + ", ".join(f"{k}={v}" for k, v in sorted(n.items())))

if __name__ == "__main__":
    main()
