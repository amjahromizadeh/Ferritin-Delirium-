library(data.table)
library(susieR)
library(coloc)
library(SNPRelate)
library(gdsfmt)
library(knitr)

## ----------------------- Inputs -----------------------
fer_file  <- "Ferritin_AF0p005.mr_ready.tsv.gz"  # <-- set to your ferritin GWAS (GRCh38)
del_file  <- "Delirium_AF0p005.mr_ready.tsv.gz"  # delirium GWAS (GRCh38)
vcf_panel <- "chr19.EUR.GRCh38.vcf.gz"           # 1KG/NYGC EUR reference panel (GRCh38)

## Region (GRCh38): APOE ±500 kb
chr_apoe   <- 19L
start_apoe <- 44421094L
end_apoe   <- 45421094L

## Sample sizes (set exact values)
N_fer <- 246000  # TODO: replace with your exact ferritin GWAS N
N_del <- 8461 + 449979

## ---------------- 1) Load & slice GWAS (both traits) ----------------
load_gwas <- function(path) {
  dt <- fread(path)
  stopifnot(all(c("SNP","chr","pos","beta","se","eaf","effect_allele","other_allele") %in% names(dt)))
  dt[, chr := as.integer(sub("^chr","", as.character(chr)))]
  dt[, `:=`(pos = as.integer(pos),
            beta = as.numeric(beta),
            se   = as.numeric(se),
            eaf  = as.numeric(eaf))]
  dt[, maf := pmin(eaf, 1 - eaf)]
  dt <- dt[
    chr == chr_apoe & pos >= start_apoe & pos <= end_apoe &
      is.finite(beta) & is.finite(se) & is.finite(maf) & maf > 0 & maf < 1
  ]
  dt <- dt[nchar(effect_allele) == 1 & nchar(other_allele) == 1]
  dt[, allele_set := paste(sort(c(effect_allele, other_allele)), collapse = "")]
  dt
}

fer <- load_gwas(fer_file)
del <- load_gwas(del_file)

cat(sprintf("Ferritin window rows: %d | Delirium window rows: %d\n", nrow(fer), nrow(del)))
stopifnot(nrow(fer) >= 10, nrow(del) >= 10)

## Remove palindromics near 0.5 (strand-ambiguous) in BOTH datasets
is_pal <- function(EA, OA) {
  (EA %in% c("A","T") & OA %in% c("A","T")) |
    (EA %in% c("C","G") & OA %in% c("C","G"))
}
fer <- fer[ !is_pal(effect_allele, other_allele) | abs(eaf - 0.5) > 0.05 ]
del <- del[ !is_pal(effect_allele, other_allele) | abs(eaf - 0.5) > 0.05 ]

## ---------------- 2) Harmonize Ferritin to Delirium EA ----------------
## join by chr:pos AND unordered allele_set to guarantee same physical allele pair
## 0) Normalize alleles and build row-wise allele sets
for (DT in list(fer, del)) {
  DT[, `:=`(
    effect_allele = toupper(effect_allele),
    other_allele  = toupper(other_allele)
  )]
}
## row-wise unordered allele pair (A/C == C/A)
fer[, aset := paste0(pmin(effect_allele, other_allele),
                     pmax(effect_allele, other_allele))]
del[, aset := paste0(pmin(effect_allele, other_allele),
                     pmax(effect_allele, other_allele))]

## 1) Primary join by chr:pos + unordered allele set
aln <- merge(
  del[, .(chr, pos,
          EA_del = effect_allele, OA_del = other_allele,
          beta_del = beta, se_del = se, eaf_del = eaf,
          aset)],
  fer[, .(chr, pos,
          EA_fer = effect_allele, OA_fer = other_allele,
          beta_fer = beta, se_fer = se, eaf_fer = eaf,
          aset)],
  by = c("chr","pos","aset")
)

cat(sprintf("Overlap by chr:pos+aset: %d SNPs\n", nrow(aln)))

## 2) Fallback: if that’s small, join by rsID (SNP), then keep DELIRIUM coords
if (nrow(aln) < 10L) {
  aln_rs <- merge(
    del[, .(SNP, chr, pos,
            EA_del = effect_allele, OA_del = other_allele,
            beta_del = beta, se_del = se, eaf_del = eaf,
            aset_del = aset)],
    fer[, .(SNP,
            EA_fer = effect_allele, OA_fer = other_allele,
            beta_fer = beta, se_fer = se, eaf_fer = eaf,
            aset_fer = aset)],
    by = "SNP"
  )
  ## keep only rows where unordered allele sets match
  aln_rs <- aln_rs[(aset_del == aset_fer)]
  
  ## build final harmonized table anchored on DELIRIUM chr/pos
  aln <- aln_rs[, .(
    chr, pos,
    EA_del, OA_del, beta_del, se_del, eaf_del,
    EA_fer, OA_fer, beta_fer, se_fer, eaf_fer,
    aset = aset_del
  )]
  cat(sprintf("Fallback by rsID matched (aset-checked): %d SNPs\n", nrow(aln)))
}

stopifnot(nrow(aln) >= 10)

## 3) Align Ferritin effect to Delirium EA
aln[, beta_fer_to_delEA :=
      fifelse(EA_fer == EA_del,  beta_fer,
              fifelse(EA_fer == OA_del, -beta_fer, NA_real_))]
aln <- aln[ is.finite(beta_del) & is.finite(se_del) &
              is.finite(beta_fer_to_delEA) & is.finite(se_fer) ]
cat(sprintf("Post-harmonization usable SNPs: %d\n", nrow(aln)))
stopifnot(nrow(aln) >= 10)
## ---------------- 3) Map to panel via chr:pos:REF:ALT ----------------
## Build panel table in the window, with an allele_set to match unordered pairs
## Build panel table in the window
gds_file <- paste0(vcf_panel, ".gds")

if (exists("genofile")) { try(closefn.gds(genofile), silent = TRUE); rm(genofile) }
if (!file.exists(gds_file)) {
  tmp <- tempfile(fileext = ".gds")
  snpgdsVCF2GDS(vcf.fn = vcf_panel, out.fn = tmp, method = "biallelic.only")
  file.rename(tmp, gds_file)
}
genofile <- openfn.gds(gds_file, readonly = TRUE, allow.duplicate = TRUE)

chr_v <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
pos_v <- read.gdsn(index.gdsn(genofile, "snp.position"))
sid_v <- read.gdsn(index.gdsn(genofile, "snp.id"))
al_v  <- read.gdsn(index.gdsn(genofile, "snp.allele"))  # "REF/ALT"

panel <- data.table(
  snp_id = as.integer(sid_v),
  chr    = as.integer(chr_v),
  pos    = as.integer(pos_v),
  allele = al_v
)[chr == chr_apoe & pos >= start_apoe & pos <= end_apoe]

alle <- tstrsplit(panel$allele, "/", fixed = TRUE)
panel[, `:=`(ref = toupper(alle[[1]]), alt = toupper(alle[[2]]))]
panel <- panel[!is.na(ref) & !is.na(alt)]

## ✅ ROW-WISE unordered allele set (this was the bug)
panel[, aset := paste0(pmin(ref, alt), pmax(ref, alt))]

## quick sanity check 
pos_only <- merge(aln[, .(chr,pos)], panel[, .(chr,pos)], by = c("chr","pos"))
cat(sprintf("QC: pos-only overlap (ignoring alleles): %d\n", nrow(pos_only)))

## Join harmonized GWAS overlap to the panel by chr:pos + aset
aln2 <- merge(
  aln[, .(chr, pos, EA_del, OA_del, beta_del, se_del, eaf_del,
          beta_fer_to_delEA, se_fer, aset)],
  panel[, .(chr, pos, ref, alt, snp_id, aset)],
  by = c("chr","pos","aset")
)

cat(sprintf("Mapped to panel (matching allele sets): %d SNPs\n", nrow(aln2)))

## Fallback: if still small, try a strand-complement aset
if (nrow(aln2) < 10L) {
  comp <- function(a) chartr("ACGT","TGCA", a)
  aln[, aset_comp := paste0(pmin(comp(EA_del), comp(OA_del)),
                            pmax(comp(EA_del), comp(OA_del)))]
  aln2b <- merge(
    aln[, .(chr, pos, EA_del, OA_del, beta_del, se_del, eaf_del,
            beta_fer_to_delEA, se_fer, aset_comp)],
    panel[, .(chr, pos, ref, alt, snp_id, aset)],
    by.x = c("chr","pos","aset_comp"),
    by.y = c("chr","pos","aset"),
    all = FALSE
  )
  if (nrow(aln2b) > nrow(aln2)) {
    aln2 <- copy(aln2b)
    cat(sprintf("Mapped to panel via strand-complement aset: %d SNPs\n", nrow(aln2)))
  }
}

stopifnot(nrow(aln2) >= 10)
## Join the already-harmonized set to the panel by chr:pos AND allele_set (unordered)
aln2 <- merge(
  aln[, .(chr, pos, EA_del, OA_del, beta_del, se_del, eaf_del,
          beta_fer_to_delEA, se_fer, aset)],
  panel[, .(chr, pos, ref, alt, snp_id, aset)],
  by = c("chr","pos","aset")
)

cat(sprintf("Mapped to panel (matching allele sets): %d SNPs\n", nrow(aln2)))
stopifnot(nrow(aln2) >= 10)

## ---------------- 4) LD on EXACT snp_id set + QC ----------------
snp_id <- aln2$snp_id

ld_obj <- snpgdsLDMat(genofile, snp.id = snp_id, method = "corr", slide = -1)
R <- ld_obj$LD

fr <- snpgdsSNPRateFreq(genofile, snp.id = snp_id)
keep <- is.finite(fr$MinorFreq) & (fr$MinorFreq > 0) & (fr$MissingRate <= 0.05)

if (sum(keep) < length(keep)) {
  R     <- R[keep, keep, drop = FALSE]
  aln2  <- aln2[keep]
  snp_id<- snp_id[keep]
}

good <- rowSums(is.finite(R)) == ncol(R) & colSums(is.finite(R)) == nrow(R)
if (sum(good) < length(good)) {
  R     <- R[good, good, drop = FALSE]
  aln2  <- aln2[good]
  snp_id<- snp_id[good]
}
stopifnot(nrow(R) >= 10 && ncol(R) == nrow(R))
R <- (R + t(R)) / 2
diag(R) <- 1
closefn.gds(genofile)

cat(sprintf("LD matrix (post-QC): %dx%d\n", nrow(R), ncol(R)))

## ---------------- 5) Z-scores & SuSiE fits ----------------
## Keep at most K SNPs by max(|z| across traits) — preserves likely signals
K <- 2000L
p <- ncol(R)
ord <- order(pmax(abs(z_del), abs(z_fer)), decreasing = TRUE)
keep <- ord[seq_len(min(p, K))]

R      <- R[keep, keep, drop = FALSE]
z_del  <- z_del[keep]
z_fer  <- z_fer[keep]
aln2   <- aln2[keep]         

## Re-sanitize R 
R <- (R + t(R)) / 2
diag(R) <- 1

cat(sprintf("Thinned to %d SNPs for SuSiE\n", ncol(R)))
Lmax_fast   <- min(5L, max(1L, ncol(R) %/% 400L))  # 3–5 is usually plenty
max_iter_fx <- 400L


if (!requireNamespace("Rfast", quietly = TRUE)) {
  message("Hint: install.packages('Rfast') to speed up large R runs.")
}

fit_del <- susie_rss(
  z = z_del, R = R, n = N_del, L = Lmax_fast,
  estimate_residual_variance = FALSE,
  max_iter = max_iter_fx,
  refine   = FALSE
)

fit_fer <- susie_rss(
  z = z_fer, R = R, n = N_fer, L = Lmax_fast,
  estimate_residual_variance = FALSE,
  max_iter = max_iter_fx,
  refine   = FALSE
)
## ---------------- 6) coloc.susie (preferred) + fallback ----------------
cs <- NULL
cs_try <- try(coloc::coloc.susie(fit_del, fit_fer), silent = TRUE)
if (!inherits(cs_try, "try-error")) cs <- cs_try

if (!is.null(cs)) {
  cat("\n=== Ferritin (GWAS) × Delirium (GWAS) — SuSiE–coloc ===\n")
  print(cs$summary)
  if (!is.null(cs$results)) {
    cat("\nTop CS-pair SNPs by SNP.PP.H4:\n")
    print(head(cs$results[order(-cs$results$SNP.PP.H4)], 10))
  } else {
    cat("\n(No per-SNP PP.H4 table returned.)\n")
  }
} else {
  message("coloc.susie failed; falling back to SuSiE×SuSiE posterior-overlap.")
  
  a1 <- fit_del$alpha; a2 <- fit_fer$alpha
  L1 <- nrow(a1);      L2 <- nrow(a2)
  stopifnot(L1 > 0, L2 > 0)
  
  ## Pairwise “shared” score = dot product of component posteriors
  PP <- a1 %*% t(a2)
  
  ## Best SNP per pair (argmax of elementwise product)
  best_idx <- vapply(seq_len(L1*L2), function(ix) {
    j <- ((ix - 1L) %% L1) + 1L
    k <- ((ix - 1L) %/% L1) + 1L
    which.max(a1[j,] * a2[k,])
  }, integer(1L))
  best_lab <- paste0(aln2$chr, ":", aln2$pos, ":", aln2$ref, ":", aln2$alt)[best_idx]
  
  top_pairs <- data.table::CJ(j = seq_len(L1), k = seq_len(L2))
  top_pairs[, `:=`(PP_shared = as.vector(PP), best_snp = best_lab)]
  
  ## Direction using most-probable SNP per component
  eg_j <- apply(a1, 1, which.max)
  eg_k <- apply(a2, 1, which.max)
  dir  <- sign(z_del[eg_j])[top_pairs$j] * sign(z_fer[eg_k])[top_pairs$k]
  top_pairs[, direction := ifelse(dir > 0, "concordant", "discordant")]
  data.table::setorder(top_pairs, -PP_shared)
  
  ## Component summaries
  comp_del <- data.table(
    component = seq_len(L1),
    lead_idx  = apply(a1, 1, which.max),
    pip_lead  = apply(a1, 1, max)
  )
  comp_del[, lead_snp := paste0(aln2$chr[lead_idx], ":", aln2$pos[lead_idx], ":", aln2$ref[lead_idx], ":", aln2$alt[lead_idx])]
  comp_del[, cs_size := sapply(susie_get_cs(fit_del, Xcorr = R)$cs, length)[component]]
  comp_del[, lead_idx := NULL]
  
  comp_fer <- data.table(
    component = seq_len(L2),
    lead_idx  = apply(a2, 1, which.max),
    pip_lead  = apply(a2, 1, max)
  )
  comp_fer[, lead_snp := paste0(aln2$chr[lead_idx], ":", aln2$pos[lead_idx], ":", aln2$ref[lead_idx], ":", aln2$alt[lead_idx])]
  comp_fer[, cs_size := sapply(susie_get_cs(fit_fer, Xcorr = R)$cs, length)[component]]
  comp_fer[, lead_idx := NULL]
  
  cat("\n### Delirium components (SuSiE)\n")
  print(kable(comp_del, format = "markdown",
              caption = "Delirium components (lead SNP, PIP, CS size)"))
  cat("\n### Ferritin components (SuSiE)\n")
  print(kable(comp_fer, format = "markdown",
              caption = "Ferritin components (lead SNP, PIP, CS size)"))
  cat("\n### Top SuSiE×SuSiE shared component pairs\n")
  print(kable(head(top_pairs, 15), format = "markdown",
              caption = "Top SuSiE×SuSiE shared component pairs"))
}