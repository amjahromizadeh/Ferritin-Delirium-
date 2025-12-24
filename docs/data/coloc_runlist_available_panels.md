# COLOC run list — **available panels only** (trimmed)

**Build convention:** Define windows in **GRCh38**; when a QTL panel is GRCh37, align by **rsID** within the same biological region.  
**Default windows:** ±250 kb (iron genes) and **±500 kb** for the **APOE** cluster (dense LD).  
**Priority:** Tier 1 > Tier 3. When two panels are listed for a gene, run both (separate COLOCs).

---

## A) Ferritin ↔ QTL

| Gene         | Tier             | QTL panel (primary)                  | QTL type | Build  | Recommended region (GRCh38)  | Window | Proxy / secondary panel(s)                                                       | Notes                               |
| ------------ | ---------------- | ------------------------------------ | -------- | ------ | ---------------------------- | ------ | -------------------------------------------------------------------------------- | ----------------------------------- |
| **APOE**     | Tier 1 (strict)  | OpenGWAS **prot-a-131**              | pQTL     | GRCh37 | chr19:44,409,976–45,409,976  | 500kb  | GTEx v10 **sQTL** Brain (e.g., Cortex) ; GTEx v8 **eQTL** Brain                  | APOE cluster (APOE/TOMM40/CEACAM19) |
| **SLC11A2**  | Tier 1 (strict)  | **eQTLGen** Whole Blood              | eQTL     | GRCh37 | chr12:50,725,057–51,225,057  | 250kb  | GTEx v8 **eQTL** Whole_Blood                                                     | Iron transporter (DMT1)             |
| **TF**       | Tier 1 (strict)  | GTEx v10 **sQTL** Liver              | sQTL     | GRCh38 | chr3:133,515,185–134,015,185 | 250kb  | GTEx v8 **eQTL** Liver ; OpenGWAS **prot-a-2950** (Serotransferrin pQTL, GRCh37) | Transferrin                         |
| **TOMM40**   | Tier 3 (lenient) | GTEx v10 **sQTL** Brain (Cortex/BA9) | sQTL     | GRCh38 | chr19:44,401,715–45,401,715  | 500kb  | GTEx v8 **eQTL** Brain                                                           | In APOE LD block                    |
| **CEACAM19** | Tier 3 (lenient) | GTEx v10 **sQTL** Brain (Cortex/BA9) | sQTL     | GRCh38 | chr19:44,119,993–45,119,993  | 500kb  | GTEx v8 **eQTL** Brain                                                           | In APOE LD block                    |
### -> Ferritin-pQTL 
#### 🧬For APOE Gene🧬
Only the primary panel used. 

**COLOC summary (APOE window, GRCh38 chr19:44.41–45.41 Mb)**
- Overlapping SNPs tested: 3,679
- Priors: p1=1e-4, p2=1e-4, p12=1e-5
- Posteriors:
    - PP.H0 (no assoc either): 5.45e-53
        
    - PP.H1 (GWAS only): 7.19e-51
        
    - PP.H2 (pQTL only): 1.55e-05
        
    - **PP.H3 (two distinct causal variants): 1.05e-03 (0.1%)**
        
    - **PP.H4 (one shared causal variant): 0.999**
        
- **Conclusion:** strong colocalization (shared variant).

PP.H4 = 0.999 → Overwhelming evidence for a single shared causal variant affecting both ferritin and APOE protein within the APOE window you defined.
PP.H3 = 0.001 → Very little support for “both traits have signals here but from different variants.”
PP.H0/H1/H2 ≈ 0 → It’s not a “no signal” region; it’s not “GWAS-only” or “pQTL-only.”
In plain words: in this 19q13 APOE region, the ferritin signal and the APOE pQTL signal are almost certainly driven by (approximately) the same genetic variant.

**Lead H4 SNPs**
- rs1065853 (PP.H4 ≈ 0.542)
- rs7412 (PP.H4 ≈ 0.458)

**Directionality sanity check**
- Unweighted sign concordance across overlapping SNPs: **46.0%** (1691/3679).
- **PP.H4-weighted sign concordance: 0.0%** → the most probable causal SNPs show **opposite** effect directions (per GWAS effect allele) for APOE protein vs ferritin.

How GPT summarized all these results: 
“COLOC between serum ferritin and APOE pQTL (prot-a-131; APOE region, ±500 kb) strongly supported a shared causal variant (PP.H4=0.999; nsnps=3,679). Posterior mass concentrated on rs1065853 and rs7412. Directionality sanity check indicated **inverse** effects of the shared allele on APOE protein vs ferritin (PP.H4-weighted sign concordance = 0%).”
### -> Ferritin-eQTL 
#### 🧬For SLC11A2 Gene🧬
Only the primary panel used. 
**region-level summary** table⬇️: 

| Locus   | Panel                 | nsnps | PP.H0 | PP.H1 | PP.H2 | PP.H3 | PP.H4 |
| :------ | :-------------------- | ----: | :---- | :---- | :---- | :---- | :---- |
| SLC11A2 | eQTLGen (Whole blood) |  1378 | 0.0%  | 0.0%  | 0.0%  | 14.7% | 85.3% |
**PP.H4 = 0.853 (~85.3%)**  
Posterior probability there is **one shared causal variant** that affects both ferritin (GWAS signal) and **SLC11A2 expression** in whole blood.  
➜ This is **strong evidence of colocalization**.

**Top 10 SNPs by SNP.PP.H4⬇️:** 

| SNP       | SNP.PP.H4 |
| :-------- | :-------- |
| rs1095998 | 18.10%    |
| rs161044  | 16.31%    |
| rs2630367 | 8.70%     |
| rs319930  | 8.46%     |
| rs150909  | 5.92%     |
| rs2046498 | 4.62%     |
| rs224564  | 3.42%     |
| rs829021  | 3.16%     |
| rs224578  | 3.00%     |
| rs706803  | 2.88%     |
Lead SNP: `rs1095998` with **SNP.PP.H4 = 0.181 (18.1%)**.  
That means, _conditional on_ H4 being true, this SNP alone carries ~18% of the “shared-causal” probability.
Top 10 SNPs (by SNP.PP.H4) together account for **~74.6%** of the shared-causal mass 
➜ The weight is **spread across several LD-correlated variants**; no single SNP dominates. This is typical when LD is moderate/high and the signals are broad.

**95% credible set size:** 21⬇️:

| SNP       | SNP.PP.H4 | Cum.PP.H4 |
| :-------- | :-------- | :-------- |
| rs1095998 | 18.10%    | 18.10%    |
| rs161044  | 16.31%    | 34.41%    |
| rs2630367 | 8.70%     | 43.11%    |
| rs319930  | 8.46%     | 51.57%    |
| rs150909  | 5.92%     | 57.49%    |
| rs2046498 | 4.62%     | 62.11%    |
| rs224564  | 3.42%     | 65.53%    |
| rs829021  | 3.16%     | 68.69%    |
| rs224578  | 3.00%     | 71.69%    |
| rs706803  | 2.88%     | 74.57%    |
| rs224566  | 2.85%     | 77.42%    |
| rs2630353 | 2.84%     | 80.25%    |
| rs224569  | 2.25%     | 82.51%    |
| rs706804  | 2.09%     | 84.60%    |
| rs2630366 | 1.95%     | 86.56%    |
| rs161045  | 1.82%     | 88.37%    |
| rs161047  | 1.49%     | 89.87%    |
| rs224570  | 1.44%     | 91.31%    |
| rs161046  | 1.37%     | 92.68%    |
| rs319938  | 1.07%     | 93.75%    |
| rs2630358 | 1.02%     | 94.77%    |
How GPT summarized all these results: 
"In the chr12 SLC11A2 window (1,378 overlapping SNPs), coloc supported a shared causal variant between ferritin GWAS and SLC11A2 expression in whole blood (eQTLGen): PP.H4 = 0.853, PP.H3 = 0.147. The shared-variant probability was distributed across multiple SNPs; the lead SNP was rs1095998 (SNP.PP.H4 = 0.181), and the top ten SNPs accounted for ~75% of the H4 mass. This indicates strong colocalization with some residual ambiguity consistent with multiple signals or LD."
##### Directionality check 🧭
What is the goal of directionality check?
After a strong COLOC (high PP.H4), we should verify that the **same allele** (e.g., the “C” allele of a SNP) pushes both traits in a **consistent direction**. For example, does the “C” allele that **lowers ferritin** also **lower SLC11A2 expression** (or raise it)? 

**How we do it:**
1. Take the **exact SNPs** used by COLOC.
2. **Align alleles** so both datasets refer to the **GWAS effect allele** (flip the eQTL sign if the counted allele differs).
3. Compare signs of the aligned eQTL effects vs GWAS betas: same sign = “agree”, opposite sign = “disagree”.
4. Summarize agreement **unweighted** (simple % of SNPs) and **weighted by SNP.PP.H4** (give more weight to SNPs that COLOC says are most likely causal).

**Our results:** 
Unweighted agreement: 
`Sign concordance (count): 833/1378 (60.4%)`

This means about **60%** of all overlapping SNPs had the **same sign** after proper allele alignment. 
Why it’s not closer to 100%: most SNPs are **not** the causal one(s); LD, noise, and cross-cohort differences can jumble per-SNP signs. This number will often be only modestly above 50%.

PP.H4-weighted agreement:
`Sign concordance (PP.H4-weighted): 100.0%` 

When we focus on the **SNPs that carry the COLOC mass** (those COLOC thinks are driving the **shared** signal), the directionality is **perfectly consistent**.
This is exactly what we want to see: the **putative shared causal signal** points both traits the **same way**.

**Lead SNP (by per-SNP H4) details:**
SNP         : rs1095998
Ferritin    : effect_allele = C, other_allele = G, eaf = 0.9543
eQTLGen     : EA (effect/assessed allele) = G, OA = C
GWAS beta   : -0.0734143   (per C allele on ferritin)
eQTL raw    : +19.8946     (this looks like a Z-score, not a beta)
Alignment   : GWAS EA = eQTL OA → we FLIP the eQTL sign
Aligned eQTL: -19.8946     (per C allele on expression) 

### Interpreting the lead SNP

- GWAS says **C** allele is associated with **lower ferritin** (beta = **–0.073**).
- After alignment, eQTL says **C** allele is associated with **lower SLC11A2 expression** (aligned slope ≈ **–19.9**, sign taken from Z-score).
    - **Same sign (both negative)** ⇒ The **C** allele moves ferritin and SLC11A2 **in the same direction** (down).
- Because **this SNP carries the largest H4 weight**, it’s a strong, coherent picture: the allele most likely to be **shared causal** reduces **both** ferritin and SLC11A2 expression. 

### -> Ferritin-sQTL  
#### 🧬For TF Gene🧬
Only the primary panel used.
Here's what the code did: 
I first subset Ferritin GWAS (GRCh38) to the TF window (chr3:133,515,185–134,015,185), computed MAF = min(EAF, 1–EAF), kept MAF ≥ 0.005, removed duplicates, and built the coloc dataset for a quantitative trait (`beta`, `se²`, `MAF`, `N`). Then, I loaded GTEx v10 Liver sQTL, identified TF via liver_genes (map gene_name == "TF" → gene_id), and kept only TF rows from liver_sig. Then, I parsed variant_id (e.g., chr3_…b38) into chr/pos, restricted to the same window, kept biallelic SNPs, and set maf = min(af, 1–af) with MAF ≥ 0.005. Then, I standardised column names to slope/slope_se; built the sQTL coloc dataset as quantitative with sdY = 1 (common practice for GTEx sQTL summary stats when per-SNP N/variance scaling isn’t explicit). 

After this, I aligned the two datasets by genomic position (same chr:pos on GRCh38), retained the intersecting SNPs (178 here), and ensured the same SNP order/IDs on both sides. Finally, I Ran `coloc.abf` with standard priors (`p1 = p2 = 1e-4`, `p12 = 1e-5`) and read the posterior probabilities (PP.H0–H4). 

The results: 
nsnps     = 178
PP.H0.abf = 1.58e-33   (neither trait associated)
PP.H1.abf = 1.37e-32   (Ferritin only)
PP.H2.abf = 1.82e-02   (TF sQTL only)
PP.H3.abf = 1.57e-01   (both traits, different causal variants)
PP.H4.abf = 8.25e-01   (both traits, shared causal variant) 

Interpretation of the numbers: 
- **nsnps = 178**: number of overlapping, QC-passing SNPs in the window used by coloc.
- **PP.H4 = 0.825** → **strong evidence of colocalisation** between the Ferritin GWAS signal and the **TF Liver sQTL** signal in the chr3 window (±250 kb). In plain words: the data support that **the same variant** is driving both the ferritin association and TF splicing in liver.
- **PP.H3 = 0.156**: some residual probability that both traits are associated **but with different variants** (e.g., if the region harbors multiple signals).
- **PP.H1/PP.H2 ≈ 0**: data do **not** favor a “one trait only” scenario.
- **PP.H0 ≈ 0**: there is clear regional signal.

Now, I want to look at the top SNP’s ferritin effect vs. sQTL slope to describe the biological direction. From the `coloc.abf` results I chose the SNP with the **largest SNP posterior for H4**, the one most likely to be the shared causal variant. Then, I joined the ferritin GWAS rows and the sQTL rows **by position** (same chr:pos on GRCh38) to make sure we were talking about exactly the **same genetic changes**. After that, I aligned alleles so effect directions are comparable.  
GTEx sQTL **slope** is per **ALT** allele; the GWAS **beta** is per **effect** allele.

- If **GWAS effect allele = sQTL ALT**, I kept the sQTL sign.
- If **GWAS effect allele = sQTL REF**, I **flipped the sQTL sign**.  
    This gave us `slope_aligned` that’s **per the same counted allele** as GWAS.

Finally, I summarized directionality: 
- Printed the **lead SNP** row with both effects and the allele-matching notes.
- Calculated **sign agreement** over all overlapping SNPs.
- Re-weighted agreement by each SNP’s **H4 posterior mass** to emphasize the variants that drive the colocalisation. Result: **100%** agreement among those.

This is the lead SNP and its information⬇️: 

|SNP       | GWAS effect_allele | GWAS other_allele | REF | ALT | EAF (GWAS)| AF (sQTL)| β_GWAS| slope_sQTL| slope_aligned|allele_match_note       |freq_support_note |
|:---------|:------------------:|:-----------------:|:---:|:---:|----------:|---------:|------:|----------:|-------------:|:-----------------------|:-----------------|
|rs3811647 |         A          |         G         |  G  |  A  |      0.279|      0.32| -0.016|     -0.699|        -0.699|GWAS EA = ALT (no flip) |EAF≈ALT           |

Here's the interpretation: 
- **Ferritin GWAS:** effect allele **A** has **beta = −0.0158** → allele **A** is associated with **lower ferritin**.
- **TF liver sQTL:** after aligning alleles to match GWAS counting, **slope_aligned = −0.699** → the **same allele A** is associated with a **lower value of the TF splicing phenotype** (e.g., less inclusion/usage of that splice junction/event).

**So:** the **same allele** (A at rs3811647) **decreases ferritin** and **decreases the TF splicing event** in liver.
🧭**Directionality:** sign concordance across overlapping SNPs was **157/178 (88.2%)**, and **PP.H4-weighted concordance was 100.0%**.

How GPT summarized all these results: 
“In the TF locus (liver), Ferritin GWAS and TF sQTL coloc support a shared causal variant (PP.H4 = 0.825; nsnps = 178). The lead SNP is rs3811647: the Ferritin-decreasing allele (A) is also associated with decreased TF splicing event usage in liver (β_ferritin = −0.016; sQTL slope_aligned = −0.699). Directionality is consistent for 88% of overlapping SNPs and for 100% of the SNP-level H4 mass.”

#### 🧬For TOMM40 Gene🧬
I couldn't run coloc for this gene using the files I was given (neither in the cortex, nor in the BA9). Here's why: 
For running sQTL coloc, I'm opening these two files (for example) for each tissue (from this folder: GTEx_Analysis_v10_sQTL): 
Brain_Cortex.v10.sQTLs.signif_pairs.parquet 
Brain_Cortex.v10.sGenes.txt.gz 

The first file contains a sparse list of only the significant results for that tissue. Each row = one variant–gene splicing pair that passed GTEx’s per-gene multiple-testing threshold. If a gene has no row here, it means no SNP passed the cutoff for that gene in this tissue (not that it wasn’t tested). 

The second file contains per-gene summary for the tissue. One row per tested gene (even if it has zero significant SNPs). This file has two important columns: pval_nominal and pval_nominal_threshold. **If pval_nominal ≤ pval_nominal_threshold ⇒ the gene should appear in signif_pairs; if min_pval_nominal > threshold ⇒ no significant pairs for that gene in this tissue.**

GTEx v10 Brain Cortex/BA9 sQTL (TOMM40): TOMM40 is present in the sGenes index (tested; gene_id = ENSG00000130204.13), but there are 0 significant variant–gene pairs in the `…sQTLs.signif_pairs.parquet` for both tissues. In the sGenes summary, the top nominal p (~2.82×10⁻⁴) is higher than the per-gene threshold (~2.90×10⁻⁵), confirming no significant sQTL at GTEx v10 power. Result: coloc could not be run using `signif_pairs` (no SNPs on the QTL side).

So, what should we do? 
We should use the all-pairs/nominal sQTL files for these tissues or just switch to brain **eQTL** for TOMM40.
I'll deal with this one later. 

**With the Proxy / secondary panel:** I used cortex eQTL. 
Running these two files from the GTEx_Analysis_v8_eQTL folder: 
Brain_Cortex.v8.signif_variant_gene_pairs.txt.gz 
Brain_Cortex.v8.egenes.txt.gz 

|          |         x|
|:---------|---------:|
|nsnps     | 1.0000000|
|PP.H0.abf | 0.8240554|
|PP.H1.abf | 0.0000035|
|PP.H2.abf | 0.1751969|
|PP.H3.abf | 0.0000000|
|PP.H4.abf | 0.0007442|
In this cortex analysis, almost all probability sits on “no GWAS signal” (H0) and “eQTL-only” (H2). The **shared signal** hypothesis (H4) is essentially **absent** (<0.1%). So there’s **no support for colocalization** of Ferritin with TOMM40 splicing eQTL in this brain tissue/window.

If you want I can coloc this gene with frontal cortex panel as well. 
#### 🧬For CEACAM19 Gene🧬
Same issue as with TOMM40 for the primary panel.
With the same cortex panel as the last gene. 
**Coloc results:**
PP.H0  ≈ 1.6e-06   (no association for either trait)
PP.H1  ≈ 3.0e-11   (Ferritin only)
PP.H2  ≈ 0.997     (eQTL only)
PP.H3  ≈ 1.5e-05   (both traits, different variants)
PP.H4  ≈ 0.00304   (both traits, same variant = colocalisation)

The posterior is **~99.7% on H2 (eQTL-only)** and **~0.3% on H4**.  
→ That means **strong evidence CEACAM19 has an eQTL in this window**, but **little/no evidence that the Ferritin GWAS shares the same causal variant** here. In short: **no colocalisation**.

If you want I can coloc this gene with frontal cortex panel as well. 

---

## B) Delirium ↔ QTL

| Gene         | Tier             | QTL panel (primary)           | QTL type | Build  | Recommended region (GRCh38)     | Window | Proxy / secondary panel(s)                       | Notes                                         |
| ------------ | ---------------- | ----------------------------- | -------- | ------ | ------------------------------- | ------ | ------------------------------------------------ | --------------------------------------------- |
| **APOE**     | Tier 1 (strict)  | OpenGWAS **prot-a-131**       | pQTL     | GRCh37 | chr19:44,421,094–45,421,094     | 500kb  | GTEx v10 **sQTL** Brain ; GTEx v8 **eQTL** Brain | APOE cluster                                  |
| **TOMM40**   | Tier 1 (strict)  | GTEx v10 **sQTL** Brain       | sQTL     | GRCh38 | chr19:44,393,408–45,393,408     | 500kb  | GoDMC **mQTL** (blood, GRCh37; proxy)            | Brain mQTL not locally available              |
| **CEACAM19** | Tier 1 (strict)  | GTEx v8 **eQTL** Brain        | eQTL     | GRCh38 | chr19:44,162,645–45,162,645     | 500kb  | GTEx v10 **sQTL** Brain                          | —                                             |
| **METTL25**  | Tier 3 (lenient) | GTEx v8 **eQTL** Brain        | eQTL     | GRCh38 | **chrNA:82,083,981–82,583,981** | 250kb  | GTEx v10 **sQTL** Brain                          | Chromosome uncertain → confirm before slicing |
| **NEBL**     | Tier 3 (lenient) | GoDMC **mQTL** (blood; proxy) | mQTL     | GRCh37 | chr10:21,019,310–21,519,310     | 250kb  | GTEx v8 **eQTL** (Brain or Whole_Blood)          | Tissue mismatch; treat as exploratory         |

### -> Delirium-pQTL
#### 🧬For APOE Gene🧬
Using the primary panel: 
- **SNPs tested:** **3,520**
- **Posterior probabilities (PP):**
    - **PP.H0 (no signals):** 6.13×10⁻¹⁹⁸ (≈0)
        
    - **PP.H1 (Delirium only):** 3.52×10⁻⁴⁸ (≈0)
        
    - **PP.H2 (APOE pQTL only):** 1.74×10⁻¹⁵⁰ (≈0)
        
    - **PP.H3 (both traits signal, different causal variants):** **~1.00**
        
    - **PP.H4 (both traits share the same causal variant):** **2.65×10⁻⁴¹ (≈0%)**

In the APOE locus, Delirium and APOE protein levels both show strong signals, but they do not share the same causal variant under this analysis. 

Since there was no evidence for colocalization with the primary pQTL panel, I'm going to test the proxy panels as well. 

Secondary panels coloc: 
1. **Delirium × GTEx v10 sQTL (Brain Cortex) for APOE**
**SNPs used:** **153** after strict, position-matched alignment/QC.

Coloc posterior probabilities (what each means)

- **PP.H0 = 1.28e−71** → Neither trait has a signal here (effectively 0%).
    
- **PP.H1 = 1.25e−50** → **Only** Delirium has a signal (≈0%).
    
- **PP.H2 = 1.03e−21** → **Only** the sQTL has a signal (≈0%).
    
- **PP.H3 = 1.00** → **Both traits have signals, but at different causal variants (≈100%).
    
- **PP.H4 = 1.74e−22** → **Both traits share the same causal variant** (≈0.000000000000000000017%).

This is a **clear “H3” outcome**: both Delirium and the Brain Cortex sQTL show association in the APOE region, but the model **strongly rejects a shared causal variant** (PP.H4 ~ 0%).

2. **Delirium × GTEx v10 sQTL (Brain BA9) for APOE***
**SNPs overlapping:** 102

- **Posteriors (coloc.abf):**
    - **PP.H0 (no signal in either):** 1.47×10⁻⁶⁴ (essentially 0)
        
    - **PP.H1 (sQTL only):** 3.77×10⁻⁴⁷ (≈0)
        
    - **PP.H2 (Delirium only):** 3.89×10⁻¹⁸ (≈0)
        
    - **PP.H3 (two distinct signals):** **~1.00**
        
    - **PP.H4 (shared causal variant):** 7.09×10⁻¹⁹ (**~0.0000000000000000007%**)

Both traits show a real association signal in this region (H0, H1, H2 are ~0), **but** the model is **overwhelmingly** confident those signals **do not come from the same causal variant** in BA9 (PP.H3 ≈ 1, PP.H4 ~ 0). 

3. Delirium × **GTEx v8 eQTL (Brain Cortex)** 
Failed attempt. 
- **What I did:**
    - Loaded the delirium GWAS and the GTEx v8 Brain Cortex **significant variant–gene pairs** + eGenes.
    - Sliced both to the **APOE ±500 kb (GRCh38)** window.
    - Tried to **align variants** two ways so coloc can use the _same SNPs_ on both sides:
        1. **By genomic position (chr:pos)** parsed from `variant_id` → **0 overlaps**.
        2. **By rsID** → the v8 “significant pairs” file **doesn’t carry any rsID column**, so rsID merging isn’t possible.
- **Why it didn’t work:**  
    Coloc needs the **exact same variants** in both datasets. The GTEx v8 file you loaded is a **sparse “significant pairs”** table (only a small set of sentinel/lead hits), and it **lacks rsIDs**. Our GWAS doesn’t necessarily include those exact chr:pos sites, so we had **no intersecting SNPs** to run coloc on.

واسه یک ژن چُسُکی 4 تا coloc کردم. همه بی‌حاصل. اَه. بیشتر از یک ساعت و نیم این چهار تا طول کشید. آخری پدرم رو درآورد. واقعا اَه. دیگه فقط پنل اول رو می‌کنم. جواب داد، داد. نداد، به کیرم. 
اصن می‌دونی چقدر برای هر کدوم از اینا برق مصرف می‌شه؟ چقدر برای هر یک coloc جی پی تی کره‌ی زمین رو گرم می‌کنه؟ الکیه مگه؟ 

Since this is a tier1 gene, I'm going to try another panel. 
Trying eQTLGen Whole Blood. 

### -> Delirium-sQTL
#### 🧬For TOMM40 Gene🧬
In GTEx v10 Brain Cortex, TOMM40 has no significant sQTLs recorded in the …sQTLs.signif_pairs.parquet file. So I'm moving on to the proxy panel(GoDMC **mQTL** (blood, GRCh37; proxy)). The mQTL panel is such a large file. It takes forever to process it. 
After very long minutes of running code and processing, guess what?
No overlapping rsIDs between GWAS (TOMM40 window, b38) and GoDMC mQTL (b37). If you want you can choose another panel for this gene and I'll run coloc.
### -> Delirium-eQTL
#### 🧬For CEACAM19 Gene🧬
Only the primary panel used. Finally some fucking results. 
**Overlap:** 5 SNPs passed QC and overlapped by position.
- **Region summary:**
    - PP.H4 (shared causal variant): **99.83%** → **strong colocalization**.
    - PP.H3 (two distinct causal variants): **0.10%**
    - PP.H0–H2: essentially ~0.
- **Lead signal:** two SNPs carry essentially all of the H4 mass:
    - **19:44646342** (PP.H4 ≈ 0.744)
    - **19:44662645** (PP.H4 ≈ 0.256)
- **95% credible set:** very likely **2 SNPs** (the two above).

Delirium × CEACAM19 (Brain Cortex v8) — coloc summary ⬇️

| nsnps| PP.H0.abf| PP.H1.abf| PP.H2.abf| PP.H3.abf| PP.H4.abf|
|-----:|---------:|---------:|---------:|---------:|---------:|
|     5|         0|         0| 0.0007387| 0.0010049| 0.9982564|
Top SNPs by PP.H4 (Delirium × CEACAM19, Cortex v8) ⬇️

|snp         |     PP.H4|
|:-----------|---------:|
|19:44646342 | 0.7436349|
|19:44662645 | 0.2563651|
|19:44677956 | 0.0000000|
|19:44683693 | 0.0000000|
|19:44619993 | 0.0000000|

95% credible set (size and cumulative PP)⬇️

| credible_set_size| cumulative_PP|
|-----------------:|-------------:|
|                 2|             1|

Lead SNP details (aligned by chr:pos key)⬇️

|SNP         | beta_gwas|   se_gwas| maf_gwas| slope_eqtl|   se_eqtl| maf_eqtl|     PP.H4|
|:-----------|---------:|---------:|--------:|----------:|---------:|--------:|---------:|
|19:44646342 | -0.116987| 0.0236655|  0.11856|  -0.585555| 0.0773764| 0.097561| 0.7436349|

##### Directionality check 🧭

Lead SNP directionality — Delirium × CEACAM19 (Cortex v8) ⬇️

|SNP         |effect_allele |other_allele |ref |alt |     eaf| af_alt| beta_gwas|     slope| slope_aligned|allele_match            |freq_support |     PP.H4|
|:-----------|:-------------|:------------|:---|:---|-------:|------:|---------:|---------:|-------------:|:-----------------------|:------------|---------:|
|19:44646342 |C             |T            |T   |C   | 0.88144|     NA| -0.116987| -0.585555|     -0.585555|GWAS EA = ALT (no flip) |NA           | 0.7436349|

1. **“Sign concordance (count): 5/5 (100.0%)”**

- We had 5 overlapping SNPs in the coloc dataset (matches our `nsnps = 5`).
- After aligning the eQTL effect to the **GWAS effect allele**, **all 5 SNPs** show the **same sign** for the GWAS effect (`beta_gwas`) and the aligned eQTL effect (`slope_aligned`).
- So, the allele that pushes the GWAS signal in one direction (e.g., increases delirium log-odds if β>0) also pushes CEACAM19 expression **in the same direction** across all overlapping SNPs. No allele-coding mishaps; internally consistent.

1. **“Sign concordance (PP.H4-weighted): 100.0%”**

- When we weight SNPs by their **per-SNP colocalization probability** (SNP.PP.H4), _all_ the H4 mass sits on SNPs whose directions agree.
- So, the SNP(s) carrying essentially all the shared-causal probability show perfectly consistent direction between the two traits.

When we look at the lead SNP table, If `beta_gwas > 0` **and** `slope_aligned > 0`: the **delirium–risk** allele (effect allele) is associated with **higher CEACAM19 expression**. 
#### 🧬For METTL25 Gene🧬
No METTL25 signif pairs in Brain Cortex v8.

Trying the secondary panel: GTEx v10 **sQTL** Brain 
No METTL25 significant sQTL pairs in Brain Cortex v10.
No METTL25 significant sQTL pairs in BA9 v10 'signif_pairs'.

Since neither the primary panel, nor the proxy panels showed any hits with this gene, I'm checking it with brain eQTLs. Starting with GTEx v8 eQTL — Brain Frontal Cortex (BA9). 
No METTL25 significant eQTL pairs in BA9 v8 'signif_variant_gene_pairs' 

Now, trying GTEx v8 eQTL — Brain Hippocampus 
No METTL25 significant eQTL pairs in Hippocampus v8 'signif_variant_gene_pairs'.

Trying Cortex v8 eQTL 
No METTL25 significant eQTL pairs in Cortex v8 'signif_variant_gene_pairs'. 

ناموسا خیلی تلاش کردم. نداره دیگه. چیکار کنه؟ 
حالا یک تیر آخری هم تو تاریکی می‌زنم. ولی واقعا آخریشه. 

Testing pQTL panel(**prot-a-131**): I don't even know if it's relevant. 
- **PP.H0 = 92.1%** → neither trait shows a (detectable) association in this window.
- **PP.H1 = 4.00%** → association for **pQTL only** (protein level), not delirium.
- **PP.H2 = 3.54%** → association for **delirium only**, not the pQTL.
- **PP.H3 = 0.154%** → both traits associate in the region but with **different causal variants**.
- **PP.H4 = 0.229%** → both associate and **share one causal variant** (this is the evidence for colocalization).

There’s **no meaningful colocalization** between METTL25 protein levels (this SomaLogic assay) and delirium risk in the tested window: **PP.H4 ≈ 0.23%** is effectively zero. Most of the probability (92%) says neither trait has a strong signal here under the model/priors with the available overlap (nsnps = **887**). 

At least I got something.
### -> Delirium-mQTL
#### 🧬For NEBL Gene🧬
Checking the primary panel: 
This mQTL file is so heavy I have to restart my laptop every time I'm done with it. 
And guess what again? 
No overlapping rsIDs between NEBL-window GWAS and GoDMC mQTL.

Now, trying with GTEx v8 Brain Cortex eQTL. 
No overlapping positions. 

Trying GTEx v8 Whole Blood eQTL. 
- **Overlap used:** `nsnps = 131` common variants in the NEBL ±250 kb window.
- **Region-level posterior:**
    - **PP.H4 = 1.39%** → weak evidence for a _shared causal variant_ affecting both delirium and NEBL expression.
    - **PP.H3 = 1.02%** → very small chance both traits have signals here but from _different_ causal variants.
    - **PP.H2 = 97.59%** → overwhelming support that **only the eQTL** (NEBL expression) shows a true association in this region; the delirium GWAS signal is essentially absent here.
    - PP.H0 and PP.H1 are ~0, meaning there is almost certainly some signal in at least one trait (it’s the eQTL).

The analysis **does not support colocalization** between delirium risk and **NEBL** whole-blood expression in this region. The data say the strong signal here is the **eQTL alone** (PP.H2 ≈ 97.6%), with little to no corresponding GWAS signal for delirium.

--------------------------------------------------------------------------
## C) Ferritin ↔ Delirium (GWAS↔GWAS) — focal regions

| Locus / Gene                 | APOE region? | Iron-axis? | Recommended region (GRCh38)  | Window |
| ---------------------------- | ------------ | ---------- | ---------------------------- | ------ |
| **APOE / TOMM40 / CEACAM19** | Yes          | No         | chr19:44,415,533–45,415,533  | 500kb  |
| **SLC11A2**                  | No           | Yes        | chr12:50,760,339–51,260,339  | 250kb  |
| **TF**                       | No           | Yes        | chr3:133,515,185–134,015,185 | 250kb  |
== APOE_TOMM40_CEACAM19: chr19:44,415,533-45,415,533 ==
PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
6.69e-260 1.14e-149 5.89e-111  1.00e+00 3.69e-112 
"PP abf for shared variant: 3.69e-110%"

== SLC11A2: chr12:50,760,339-51,260,339 ==
PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
 4.69e-19  5.86e-01  3.06e-19  3.82e-01  3.19e-02 
"PP abf for shared variant: 3.19%"

== TF: chr3:133,515,185-134,015,185 ==
PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
0.1150    0.5620    0.0517    0.2530    0.0190 
1 "PP abf for shared variant: 1.9%"

Results by region:
1) **APOE/TOMM40/CEACAM19 (chr19:44.42–45.42 Mb; 500 kb)**
`PP.H0 ~ 6.7e-260 | PP.H1 ~ 1.1e-149 | PP.H2 ~ 5.9e-111 | PP.H3 ≈ 1.000 | PP.H4 ~ 3.7e-112`
Interpretation: Extremely strong evidence that both traits have regional signals, but not the same variant (H3 ~ 1.00, H4 ≈ 0).
Conclusion: No colocalisation here; ferritin and delirium signals are distinct in the APOE region.

2) **SLC11A2 (chr12:50.76–51.26 Mb; 250 kb)**
`PP.H0 ~ 4.7e-19 | PP.H1 = 0.586 | PP.H2 ~ 3.1e-19 | PP.H3 = 0.382 | PP.H4 = 0.032`
Interpretation: The most likely scenario is ferritin-only (H1 ≈ 59%). There’s some chance both traits are associated but with different variants (H3 ≈ 38%). Shared-variant (H4) is low (3.2%).
Conclusion: No evidence of colocalisation. At this locus, ferritin likely has a signal; delirium does not (or, if it does, it’s not the same variant).

3) **TF (chr3:133.52–134.02 Mb; 250 kb)**
`PP.H0 = 0.115 | PP.H1 = 0.562 | PP.H2 = 0.052 | PP.H3 = 0.253 | PP.H4 = 0.019`
Interpretation: Again points to ferritin-only (H1 ≈ 56%). Some probability of both traits but different variants (H3 ≈ 25%), and very low shared-variant support (H4 ≈ 2%).
Conclusion: No colocalisation at TF; ferritin likely drives the regional association.

### Notes & caveats
- If a QTL panel is **GRCh37** (eQTLGen, OpenGWAS pQTL, GoDMC), align by **rsID** lists defined from your GRCh38 region, or liftover a minimal SNP list if needed.
- For APOE cluster, run **multiple brain tissues** (Cortex, BA9, Hippocampus) to stabilize COLOC. Keep both eQTL and sQTL runs if feasible.
- Document sQTL “significant-only” limitation (GTEx sQTL tar). If a window is sparse, supplement with eQTL or expand ±100–200 kb as sensitivity.
