# Ferritin–Delirium Genomic Analyses

Code and scripts for the integrative genomic study of serum ferritin and delirium risk, including  two-sample Mendelian randomization, SMR/HEIDI, colocalisation, and SuSiE-based fine-mapping. 

## 1. Study overview

This repository contains the analysis code used in the study:

**​Serum ferritin and delirium risk: Integrative genomic analysis identifies locus-specific signals at 19q13 with no genome-wide causal effect**

The main goals are to:
- test whether genetically proxied serum ferritin has a causal effect on delirium risk using two-sample Mendelian randomization (MR);
- characterise cis-mediated molecular mechanisms using SMR/HEIDI across blood and brain QTL resources;
- evaluate regional co-localisation of ferritin, delirium, and some specific regions QTL signals;
- perform SuSiE-based fine-mapping and SuSiE–coloc in selected regions (e.g. APOE locus).

All scripts are written in R and are designed to be run on GWAS/QTL summary statistics obtained from public resources (not redistributed in this repository).

## 2. Data and inputs

The analyses use the following external datasets (summary statistics):

- **Serum ferritin GWAS**  
  - Source: [GWAS catalog](https://www.ebi.ac.uk/gwas/studies/GCST90270865)   

- **Delirium GWAS**  
  - Source: [GWAS catalog](https://www.ebi.ac.uk/gwas/studies/GCST90473243)

- **QTL panels** 
  - pQTL: OpenGWAS **prot-a-131**
  - sQTL: 
	  - GTEx v10 sQTL Liver
	  - GTEx v10 sQTL Brain (Cortex/BA9)
  - eQTL: 
	  - eQTLGen Whole Blood
	  - GTEx v8 eQTL Liver
	  - GTEx v8 eQTL Brain
	  - GTEx v8 eQTL Whole_Blood
  - mQTL: GoDMC mQTL 

- **LD reference panel**  
  - 1000 Genomes EUR (PLINK and/or GDS format), used for SMR, coloc, and SuSiE–coloc.

**Note:** GWAS and QTL summary statistics are not redistributed in this repository. Users must obtain them directly from the original sources and adjust paths in the scripts accordingly.

## 3. Repository structure

The repository is organised as:

- `code/` – all analysis scripts (R).
- `outputs/` – key result tables produced by the scripts (e.g. MR summary tables, SMR/HEIDI hits, coloc/SuSiE output, analysis logs).
- `docs/` – figures, schematic diagrams, or additional documentation.
- `LICENSE` – license for this repository.
- `README.md` – this file.

### 3.1 Script mapping

A detailed mapping of scripts to their purpose is provided in the following table:

| Index | Filename                                              | The code it includes                                                                                                                                                                                                                                |
| ----- | ----------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1     | Ferritin-Delirium Two sample MR                       | Two-sample MR pipeline for Ferritin-Delirium (instrument selection, harmonisation, MR + sensitivity analyses, MR-PRESSO, and all main/supplementary MR plots).                                                                                      |
| 2     | Ferritin SMR results                                  | Aggregates ferritin SMR/HEIDI results across blood, brain, and liver panels, flags genome-wide SMR+HEIDI hits, classifies genes by tissue sharing, generates target gene dossiers, and prepares a −log₁₀(p_SMR) heatmap with HEIDI failures marked. |
| 3     | Bonferroni Threshold on SMR results                   | Applies Bonferroni correction to SMR results, flags significant genes, and writes filtered SMR/HEIDI tables.                                                                                                                                        |
| 4     | Ferritin Delirium coloc                               | Performs regional colocalisation (coloc.abf) between serum ferritin and delirium GWAS around 3 predefined loci and outputs PP.H0–H4 summaries per gene window.                                                                                      |
| 5     | Delirium × GTEx v10 sQTL (Brain Cortex) (APOE region) | Colocalisation of Delirium GWAS with GTEx v10 Brain Cortex sQTLs in the APOE window. (method: coloc.abf)                                                                                                                                            |
| 6     | Delirium x eQTL (whole blood)(APOE region)            | Colocalisation of Delirium GWAS with eQTLGen whole-blood APOE eQTLs. (method: coloc.abf)                                                                                                                                                            |
| 7     | Delirium x eQTL(Brain Cortex)(CEACAM19 region)        | Colocalisation of Delirium GWAS with brain cortex eQTLGen in the CEACAM19 window. (method: coloc.abf)                                                                                                                                               |
| 8     | Delirium x pQTL(APOE region)                          | Colocalisation of Delirium GWAS with pQTL in the APOE  window. (method: coloc.abf)                                                                                                                                                                  |
| 9     | Delirium x pQTL(METTL25 region)                       | Colocalisation of Delirium GWAS with pQTL in the METTL25 window. (method: coloc.abf)                                                                                                                                                                |
| 10    | Ferritin x eQTL(Brain Cortex)(TOMM40 region)          | Colocalisation of Ferritin GWAS with brain cortex eQTLGen the TOMM40 window. (method: coloc.abf)                                                                                                                                                    |
| 11    | Ferritin x eQTL(whole blood)(SLC11A2 region)          | Colocalisation of Ferritin GWAS with eQTLGen whole-blood SLC11A2 eQTLs. (method: coloc.abf)                                                                                                                                                         |
| 12    | Ferritin x pQTL(APOE region)                          | Colocalisation of Ferritin GWAS with pQTL in the APOE  window. (method: coloc.abf)                                                                                                                                                                  |
| 13    | Ferritin x sQTL(TF, TOMM40 & CEACAM19 regions)        | Colocalisation of Ferritin GWAS with sQTL in the TF, TOMM40, and CEACAM19 windows. (method: coloc.abf)                                                                                                                                              |
| 14    | SuSiE–coloc Delirium x APOE pQTL                      | SuSiE fine-mapping and `coloc.susie` colocalisation between Delirium GWAS and pQTL (prot-a-131) in the APOE  window.                                                                                                                                |
| 15    | SuSiE–coloc Delirium x APOE sQTL                      | SuSiE fine-mapping and `coloc.susie` colocalisation between Delirium GWAS and sQTL in the APOE  window.                                                                                                                                             |
| 16    | SuSiE-coloc Ferritin x APOE pQTL                      | SuSiE fine-mapping and `coloc.susie` colocalisation between Ferritin GWAS and pQTL (prot-a-131) in the APOE window.                                                                                                                                 |
| 17    | SuSiE-coloc Ferritin x Delirium                       | SuSiE fine-mapping and `coloc.susie` colocalisation between Ferritin and Delirium GWAS in the APOE region.                                                                                                                                          |
For a quick overview, the main **analysis blocks** are:

- **Two-sample MR:**  
  - `Ferritin-Delirium Two sample MR.R`

- **SMR/HEIDI processing and filtering:**  
  - `Ferritin SMR results.R`, `Bonferroni Threshold on SMR results.R`

- **Colocalisation and fine-mapping:**  
  - `Ferritin Delirium coloc.R`
  - `SuSiE–coloc Delirium x APOE pQTL.R`
  - `SuSiE-coloc Ferritin x APOE pQTL.R`
  - `SuSiE-coloc Ferritin x Delirium (GWAS).R`

- **Region-specific QTL analyses:**  
  - `Delirium × GTEx v10 sQTL (Brain Cortex) (APOE region).R`  
  - `Delirium x eQTL (whole blood)(APOE region).R`  
  - `Delirium x eQTL(Brain Cortex)(CEACAM19 region).R`

## 4. Software requirements

Analyses were performed in **R (4.5.1)** with the following key packages:

- `data.table`
- `TwoSampleMR`
- `MRPRESSO`
- `coloc`
- `susieR`
- `SNPRelate`
- `gdsfmt`
- `knitr` / `rmarkdown`

We recommend installing packages from CRAN and, where needed, from Bioconductor or GitHub. 

## 5. Outputs 

- **Two-sample MR script**
  - Outputs: MR estimates.md [[MR estimates]]

- **Coloc and SuSiE–coloc scripts**
  - Outputs: [[Coloc outputs]] [[Susie; process and outputs]]
## 6. Citation

If you use this code, please cite:

> Full manuscript citation here (once published).  

and, where appropriate, the original GWAS and QTL sources listed in the Data section.

## 7. License

The code in this repository is released under the **[MIT]** license.

This license applies to the code only and **does not** cover the external GWAS/QTL summary statistics, which remain under the terms specified by their original providers.


