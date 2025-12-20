Ferritin and Delirium Two-sample MR. 
Datasets used: 
Ferritin_AF0p005.mr_ready.tsv.gz 
Delirium_AF0p005.mr_ready.tsv.gz 

I did the Two-sample MR test and here are the results broken down to as much details as possible: 
1. **mean_F** = 78.76 
	For each instrument SNP we compute an **F-statistic** for its association with the **exposure**. 
	F measures **instrument strength**. Weak instruments bias MR toward the confounded observational association and inflate Type I error. **F ≥ 10** for each SNP is “strong enough” by classic criteria. - `mean_F = 78.76` means that, **on average**, our instruments are **very strong** (well above 10). Our SNP set has plenty of predictive signal for ferritin; weak-IV bias is unlikely to explain our MR results. 
--------------------------------------------------------------------------
2. **main causal estimates** by method.⬇️

| id.exposure | id.outcome | outcome  | exposure | method                                                    | nsnp |         b |        se |      pval |
| :---------- | :--------- | :------- | :------- | :-------------------------------------------------------- | ---: | --------: | --------: | --------: |
| AlghdU      | TBMBDR     | Delirium | Ferritin | Inverse variance weighted                                 |   70 | 0.0822032 | 0.0764773 | 0.2824330 |
| AlghdU      | TBMBDR     | Delirium | Ferritin | Inverse variance weighted (multiplicative random effects) |   70 | 0.0822032 | 0.0764773 | 0.2824330 |
| AlghdU      | TBMBDR     | Delirium | Ferritin | MR Egger                                                  |   70 | 0.1617036 | 0.1414843 | 0.2570848 |
| AlghdU      | TBMBDR     | Delirium | Ferritin | Weighted median                                           |   70 | 0.0898712 | 0.0974262 | 0.3562921 |

Our numbers: 
- IVW (and IVW-MRE): **b = 0.0822**, SE **0.0765**, p **0.282** → **not significant**.
- MR-Egger: **b = 0.162**, p **0.257** → also **not significant**.
- Weighted median: **b = 0.0899**, p **0.356** → **not significant**.

Across standard estimators, there’s **no robust evidence** that raising ferritin changes delirium risk. 

--------------------------------------------------------------------------
3. This table gives **heterogeneity** tests (Cochran’s **Q**) for IVW and MR-Egger. ⬇️

|id.exposure |id.outcome |outcome  |exposure |method                    |        Q| Q_df|    Q_pval|
|:-----------|:----------|:--------|:--------|:-------------------------|--------:|----:|---------:|
|AlghdU      |TBMBDR     |Delirium |Ferritin |MR Egger                  | 96.36125|   68| 0.0134358|
|AlghdU      |TBMBDR     |Delirium |Ferritin |Inverse variance weighted | 96.99547|   69| 0.0147892|
**Columns:**
- `method` – IVW or MR-Egger.
- `Q` – Cochran’s Q statistic (measures **excess dispersion** of per-SNP causal estimates).
- `Q_df` – degrees of freedom (**nSNP−1** for IVW; **nSNP−2** for Egger).
- `Q_pval` – p-value testing **“no heterogeneity”**. 

**Interpretation:**
- **Significant heterogeneity** (p≈0.015) → SNP estimates vary **more than expected** by chance.
- This can reflect **pleiotropy**, measurement error, or mixed mechanisms.
--------------------------------------------------------------------------
4. **MR-Egger intercept** test for **directional (unbalanced) horizontal pleiotropy** ⬇️

| id.exposure | id.outcome | outcome  | exposure | egger_intercept |        se |      pval |
| :---------- | :--------- | :------- | :------- | --------------: | --------: | --------: |
| AlghdU      | TBMBDR     | Delirium | Ferritin |      -0.0034819 | 0.0052046 | 0.5057626 |

**Interpretation:**
- **No evidence** of directional pleiotropy on average.
- Combined with Q>0, this suggests **heterogeneity without a consistent direction**, i.e., **balanced** or idiosyncratic pleiotropy.
--------------------------------------------------------------------------
5. steiger directionality test ⬇️ 

| id.exposure | id.outcome | exposure | outcome  | snp_r2.exposure | snp_r2.outcome | correct_causal_direction | steiger_pval |
| :---------- | :--------- | :------- | :------- | --------------: | -------------: | :----------------------- | -----------: |
| AlghdU      | TBMBDR     | Ferritin | Delirium |       0.0206426 |      0.0002185 | TRUE                     |            0 |
Steiger compares how much variance the IVs explain in **exposure** vs **outcome**.

- `snp_r2.exposure` – total R2R^2R2 (approx.) in **ferritin** explained by the IVs.
- `snp_r2.outcome` – total R2R^2R2 in **delirium** explained by the same IVs.
- `correct_causal_direction` – **TRUE** if RX2>RY2R^2_{X} > R^2_{Y}RX2​>RY2​, supporting **Ferritin → Delirium** rather than reverse.
- `steiger_pval` – p-value that this ordering is genuine. 

Our numbers:

- RX2=0.02064R^2_X = 0.02064RX2​=0.02064 (≈ **2.06%** variance in ferritin),
- RY2=0.0002185R^2_Y = 0.0002185RY2​=0.0002185 (≈ **0.022%** in delirium),
- `correct_causal_direction = TRUE`, `steiger_pval = 0`.

**Interpretation:**

- The instruments “explain” much more of ferritin than delirium—**direction is coherent** with our MR (Ferritin → Delirium).
- This does **not** prove causality by itself; it just says reverse causation is unlikely to explain our MR estimate.
--------------------------------------------------------------------------
6. MR-PRESSO ⬇️ 

|id.exposure |id.outcome |exposure |outcome  | snp_r2.exposure| snp_r2.outcome|correct_causal_direction | steiger_pval|
|:-----------|:----------|:--------|:--------|---------------:|--------------:|:------------------------|------------:|
|AlghdU      |TBMBDR     |Ferritin |Delirium |       0.0206426|      0.0002185|TRUE                     |            0|
> # print markdown table
> kable(mrp, format = "markdown")

|Exposure      |MR Analysis       | Causal Estimate|        Sd|   T-stat|   P-value|
|:-------------|:-----------------|---------------:|---------:|--------:|---------:|
|beta.exposure |Raw               |       0.0878940| 0.0762235| 1.153109| 0.2527882|
|beta.exposure |Outlier-corrected |       0.1198282| 0.0585500| 2.046597| 0.0445093|

|        x|
|--------:|
| 101.5063|

|     x|
|-----:|
| 0.017|

|    RSSobs|Pvalue |
|---------:|:------|
| 0.0000734|1      |
| 0.0000066|1      |
| 0.0000379|1      |
| 0.0000045|1      |
| 0.0006172|1      |
| 0.0000222|1      |
| 0.0014621|1      |
| 0.0000307|1      |
| 0.0000582|1      |
| 0.0003319|1      |
| 0.0003548|1      |
| 0.0003439|1      |
| 0.0000020|1      |
| 0.0000372|1      |
| 0.0003861|1      |
| 0.0000017|1      |
| 0.0016108|0.852  |
| 0.0002910|1      |
| 0.0000146|1      |
| 0.0001093|1      |
| 0.0038807|1      |
| 0.0076431|0.213  |
| 0.0009212|1      |
| 0.0001077|1      |
| 0.0001705|1      |
| 0.0000209|1      |
| 0.0002058|1      |
| 0.0000476|1      |
| 0.0001044|1      |
| 0.0024310|1      |
| 0.0000462|1      |
| 0.0001840|1      |
| 0.0001392|1      |
| 0.0006374|1      |
| 0.0000036|1      |
| 0.0000326|1      |
| 0.0000605|1      |
| 0.0000001|1      |
| 0.0004375|1      |
| 0.0003483|1      |
| 0.0000000|1      |
| 0.0002588|1      |
| 0.0001543|1      |
| 0.0001011|1      |
| 0.0005950|1      |
| 0.0000507|1      |
| 0.0003182|1      |
| 0.0001839|1      |
| 0.0008182|1      |
| 0.0000430|1      |
| 0.0000417|1      |
| 0.0001174|1      |
| 0.0000648|1      |
| 0.0002504|1      |
| 0.0001121|1      |
| 0.0001313|1      |
| 0.0010630|1      |
| 0.0024048|1      |
| 0.0387929|<0.071 |
| 0.0000307|1      |
| 0.0017541|1      |
| 0.0011842|1      |
| 0.0001738|1      |
| 0.0000170|1      |
| 0.0000104|1      |
| 0.0009799|1      |
| 0.0002285|1      |
| 0.0000007|1      |
| 0.0002711|1      |
| 0.0001086|1      |
| 0.0000136|1      |

|  x|
|--:|
| 59|

|              |         x|
|:-------------|---------:|
|beta.exposure | -26.65004|

|     x|
|-----:|
| 0.705|

Our output:

- **Global test p = 0.017** → **at least one outlier** is present.
- **Outlier test** table: lists each SNP’s outlier statistic. Only the **index = 59** is flagged downstream.
- **Main MR results**:
    
    - **Raw**: β **0.0879**, SE **0.0762**, p **0.253** (not significant).
        
    - **Outlier-corrected**: β **0.1198**, SE **0.0586**, p **0.0445** (nominally significant **after** removing outlier(s)).
        
- **Distortion test p = 0.705** → the change from raw → corrected is **not statistically significant**.  
    In words: we did find an outlier, and removing it nudged the estimate to nominal significance, but that shift isn’t itself significant; treat the “corrected” p=0.045 as fragile/exploratory.


