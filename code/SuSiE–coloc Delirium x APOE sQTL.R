library(data.table)
library(arrow)
library(susieR)
library(coloc)
library(SNPRelate)
library(gdsfmt)

## ---- file inputs ----
del_file    <- "Delirium_AF0p005.mr_ready.tsv.gz"
sqtl_file   <- "Brain_Cortex.v10.sQTLs.signif_pairs.parquet"
sgenes_file <- "Brain_Cortex.v10.sGenes.txt.gz"
vcf_file    <- "chr19.EUR.GRCh38.vcf.gz"

## ---- region (GRCh38): APOE ±500 kb ----
chr_apoe   <- 19L
start_apoe <- 44421094L
end_apoe   <- 45421094L

## ---- sample sizes ----
N_del  <- 8461 + 449979      
N_sqtl <- 838                

## ===================== 1) LOAD & SLICE DELIRIUM GWAS ======================
del <- fread(del_file)
stopifnot(all(c("SNP","chr","pos","beta","se","eaf","effect_allele","other_allele") %in% names(del)))

del[, chr := as.integer(sub("^chr","", as.character(chr)))]
del[, `:=`(pos = as.integer(pos), beta = as.numeric(beta), se = as.numeric(se), eaf = as.numeric(eaf))]
del[, maf := pmin(eaf, 1 - eaf)]

gwas <- del[
  chr == chr_apoe & pos >= start_apoe & pos <= end_apoe &
    !is.na(beta) & !is.na(se) & !is.na(maf) & maf > 0 & maf < 1
]
gwas <- unique(gwas, by = "SNP")
gwas <- gwas[nchar(effect_allele) == 1 & nchar(other_allele) == 1]
cat(sprintf("GWAS window rows: %d (unique rsIDs)\n", nrow(gwas)))

## ===================== 2) LOAD APOE sQTL (GTEx v10, Cortex) ===============
genes <- fread(sgenes_file)  # Brain_Cortex.v10.sGenes.txt.gz
stopifnot(all(c("gene_name","gene_id") %in% names(genes)))
genes[, gene_base := sub("\\.\\d+$","", gene_id)]
apoe_ids <- unique(genes[gene_name == "APOE", gene_base])
stopifnot(length(apoe_ids) > 0)

sqtl <- as.data.table(read_parquet(sqtl_file))  # Brain_Cortex.v10.sQTLs.signif_pairs.parquet
need_cols <- c("variant_id","slope","slope_se")
stopifnot(all(need_cols %in% names(sqtl)))

gene_col <- intersect(names(sqtl), c("molecular_trait_id","phenotype_id","gene_id"))[1]
stopifnot(length(gene_col) == 1)

## robust Ensembl extraction (drop version)
sqtl[, gene_base := {
  x <- get(gene_col)
  m <- regexpr("ENSG\\d+(?:\\.\\d+)?", x, perl = TRUE)
  y <- ifelse(m > 0, substr(x, m, m + attr(m, "match.length") - 1), NA_character_)
  sub("\\.\\d+$", "", y)
}]

## First try gene-based subset 
sqtl_apoe <- sqtl[gene_base %in% apoe_ids & !is.na(slope) & !is.na(slope_se)]

## If empty (shouldn’t be now), fall back to locus-based within window
if (nrow(sqtl_apoe) == 0L) {
  parts <- tstrsplit(sqtl$variant_id, "_", fixed = TRUE)
  stopifnot(length(parts) >= 2)
  sqtl[, chr := as.integer(sub("^chr","", parts[[1]]))]
  sqtl[, pos := as.integer(parts[[2]])]
  sqtl_apoe <- sqtl[
    chr == chr_apoe & pos >= start_apoe & pos <= end_apoe &
      !is.na(slope) & !is.na(slope_se)
  ]
}
stopifnot(nrow(sqtl_apoe) > 0)

## Parse variant_id -> chr/pos/ref/alt
parts <- tstrsplit(sqtl_apoe$variant_id, "_", fixed = TRUE)
stopifnot(length(parts) >= 5L)
sqtl_apoe[, `:=`(
  chr = as.integer(sub("^chr","", parts[[1]])),
  pos = as.integer(parts[[2]]),
  ref = parts[[3]],
  alt = parts[[4]]
)]
sqtl_apoe <- sqtl_apoe[
  chr == chr_apoe & pos >= start_apoe & pos <= end_apoe &
    nchar(ref) == 1 & nchar(alt) == 1
]
cat(sprintf("APOE sQTL rows (Cortex v10) in window: %d\n", nrow(sqtl_apoe)))
stopifnot(nrow(sqtl_apoe) > 0)

## ===================== 3) OPEN GDS & BUILD PANEL MAP ======================
gds_file <- paste0(vcf_file, ".gds")
if (!file.exists(gds_file)) {
  tmp <- tempfile(fileext = ".gds")
  snpgdsVCF2GDS(vcf.fn = vcf_file, out.fn = tmp, method = "biallelic.only")
  file.rename(tmp, gds_file)
}
genofile <- snpgdsOpen(gds_file)

chr_v <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
pos_v <- read.gdsn(index.gdsn(genofile, "snp.position"))
rs_v  <- read.gdsn(index.gdsn(genofile, "snp.rs.id"))
sid_v <- read.gdsn(index.gdsn(genofile, "snp.id"))
al_v  <- read.gdsn(index.gdsn(genofile, "snp.allele")) # "A/G" (REF/ALT)

panel <- data.table(
  snp_id = as.integer(sid_v),
  rsid   = rs_v,
  chr    = as.integer(chr_v),
  pos    = as.integer(pos_v),
  allele = al_v
)[chr == chr_apoe & pos >= start_apoe & pos <= end_apoe]

alle <- tstrsplit(panel$allele, "/", fixed = TRUE)
panel[, `:=`(ref = alle[[1]], alt = alle[[2]])]
panel <- panel[!is.na(ref) & !is.na(alt)]

## ===================== 4) MAP sQTL -> PANEL (handle swapped alleles) ======
## direct REF/ALT match
m1 <- merge(
  sqtl_apoe[, .(chr, pos, ref, alt, slope, slope_se)],
  panel[,  .(chr, pos, ref, alt, rsid, snp_id)],
  by = c("chr","pos","ref","alt")
)[, `:=`(slope_eff = slope)]

## swapped (panel has ref/alt swapped vs sQTL): flip sign
m2_in <- copy(sqtl_apoe)[, c("ref","alt") := .(alt, ref)]
m2 <- merge(
  m2_in[, .(chr, pos, ref, alt, slope, slope_se)],
  panel[, .(chr, pos, ref, alt, rsid, snp_id)],
  by = c("chr","pos","ref","alt")
)[, `:=`(slope_eff = -slope)]

sqtl_map <- unique(rbindlist(list(m1, m2), use.names = TRUE, fill = TRUE), by = "snp_id")
cat(sprintf("sQTL variants mapped to panel (unique snp_id): %d\n", nrow(sqtl_map)))
stopifnot(nrow(sqtl_map) > 0)

## ===================== 5) MERGE WITH GWAS BY chr:pos & HARMONISE ==========
## Join by genomic position to avoid rsID dropouts
gwas_pos <- gwas[, .(chr, pos, SNP, effect_allele, other_allele,
                     beta_gwas = beta, se_gwas = se, eaf, maf)]

aln <- merge(
  gwas_pos,
  sqtl_map[, .(chr, pos, snp_id, rsid, ref, alt, slope_eff, slope_se)],
  by = c("chr","pos")
)

## Palindromics: keep only if |EAF - 0.5| > 0.05
is_pal <- ( (aln$effect_allele %in% c("A","T") & aln$other_allele %in% c("A","T")) |
              (aln$effect_allele %in% c("C","G") & aln$other_allele %in% c("C","G")) )
aln <- aln[ !is_pal | abs(eaf - 0.5) > 0.05 ]

## Require EA/OA to match panel REF/ALT in some orientation
aln <- aln[
  (effect_allele == alt & other_allele == ref) |
    (effect_allele == ref & other_allele == alt)
]

## Align sQTL slope to GWAS EA (panel is REF/ALT; slope_eff already swap-corrected)
aln[, slope_aligned := fifelse(effect_allele == alt,  slope_eff,
                               fifelse(effect_allele == ref, -slope_eff, NA_real_))]

## Final QC
aln <- aln[ is.finite(beta_gwas) & is.finite(se_gwas) &
              is.finite(slope_aligned) & is.finite(slope_se) ]

cat(sprintf("Aligned, harmonized overlap (by chr:pos): %d SNPs\n", nrow(aln)))
stopifnot(nrow(aln) >= 10)

## ===================== 6) BUILD LD FOR EXACT snp_id SET ===================
snp_id <- aln$snp_id

ld_obj <- snpgdsLDMat(genofile, snp.id = snp_id, method = "corr", slide = -1)
R <- ld_obj$LD

## freq/callrate QC on the same SNPs
fr <- snpgdsSNPRateFreq(genofile, snp.id = snp_id)
keep <- is.finite(fr$MinorFreq) & (fr$MinorFreq > 0) & (fr$MissingRate <= 0.05)

if (sum(keep) < length(keep)) {
  R      <- R[keep, keep, drop = FALSE]
  aln    <- aln[keep]
  snp_id <- snp_id[keep]
}

## drop any residual non-finite rows/cols (minimal cut)
good <- rowSums(is.finite(R)) == ncol(R) & colSums(is.finite(R)) == nrow(R)
if (sum(good) < length(good)) {
  R      <- R[good, good, drop = FALSE]
  aln    <- aln[good]
  snp_id <- snp_id[good]
}

stopifnot(nrow(R) >= 10 && ncol(R) == nrow(R))
R <- (R + t(R)) / 2
diag(R) <- 1
cat(sprintf("LD matrix (post-QC): %dx%d\n", nrow(R), ncol(R)))

## ===================== 7) Z-SCORES ========================================
z1 <- aln$beta_gwas     / aln$se_gwas
z2 <- aln$slope_aligned / aln$slope_se
stopifnot(length(z1) == nrow(R) && length(z2) == nrow(R))

## ===================== 8) RUN SUSIE & COLOC.SUSIE =========================
Lmax <- max(5L, min(10L, ncol(R) %/% 50))

# label R by SNPs so outputs carry names
snp_names <- if (!is.null(gwas$SNP) && length(gwas$SNP) == nrow(R)) gwas$SNP else {
  if (!is.null(aln$SNP) && length(aln$SNP) == nrow(R)) aln$SNP else
    if (!is.null(colnames(R))) colnames(R) else paste0("s", seq_len(nrow(R)))
}
dimnames(R) <- list(snp_names, snp_names)

# SuSiE fits (match your first test: fixed residual variance)
fit1 <- susie_rss(z = z1, R = R, n = N_del,  L = Lmax,
                  estimate_residual_variance = FALSE, max_iter = 2000, refine = TRUE)
fit2 <- susie_rss(z = z2, R = R, n = N_sqtl, L = Lmax,
                  estimate_residual_variance = FALSE, max_iter = 2000, refine = TRUE)

# alpha: L x p posterior allocations; coerce & QC
a1 <- as.matrix(fit1$alpha)
a2 <- as.matrix(fit2$alpha)
if (is.null(dim(a1)) || is.null(dim(a2))) stop("SuSiE returned empty alpha; check fits.")
if (ncol(a1) != nrow(R) || ncol(a2) != nrow(R)) stop("alpha columns not aligned to R; check inputs.")
a1[!is.finite(a1)] <- 0
a2[!is.finite(a2)] <- 0

L1 <- nrow(a1); L2 <- nrow(a2)

# Pairwise shared-causal probability for components (j,k): sum_s a1[j,s]*a2[k,s]
# Compute all pairs at once (L1 x L2)
S <- a1 %*% t(a2)   # rows: trait1 components, cols: trait2 components
S[!is.finite(S)] <- 0

# Best SNP per pair = argmax over s of a1[j,s] * a2[k,s]
best_idx <- matrix(NA_integer_, nrow = L1, ncol = L2)
for (j in seq_len(L1)) {
  aj <- a1[j, ]
  for (k in seq_len(L2)) {
    bk <- a2[k, ]
    w  <- aj * bk
    best_idx[j, k] <- if (all(!is.finite(w) | w <= 0)) NA_integer_ else which.max(w)
  }
}
best_snp <- matrix(NA_character_, nrow = L1, ncol = L2)
best_snp[!is.na(best_idx)] <- snp_names[na.omit(as.vector(best_idx))]

# Credible set sizes (using Xcorr=R to be safe)
cs1 <- susie_get_cs(fit1, Xcorr = R)
cs2 <- susie_get_cs(fit2, Xcorr = R)
cs1_size <- rep(NA_integer_, L1)
cs2_size <- rep(NA_integer_, L2)
if (!is.null(cs1$cs_index)) {
  for (i in seq_along(cs1$cs_index)) cs1_size[cs1$cs_index[i]] <- length(cs1$cs[[i]])
}
if (!is.null(cs2$cs_index)) {
  for (i in seq_along(cs2$cs_index)) cs2_size[cs2$cs_index[i]] <- length(cs2$cs[[i]])
}

# Tidy table
library(data.table)
pair <- CJ(j = seq_len(L1), k = seq_len(L2))
pair[, `:=`(
  PP_shared = as.vector(S),
  best_snp  = as.vector(best_snp),
  cs1_size  = cs1_size[j],
  cs2_size  = cs2_size[k]
)]
setorder(pair, -PP_shared)

cat("\n=== SuSiE × SuSiE (fallback) — top component pairs by sharing prob ===\n")
print(pair[1:min(15, .N)])

cat(sprintf("\nTrait1 (GWAS): %d components; Trait2 (sQTL): %d components\n", L1, L2))



# --- 1) Collect helpers from current fits ---
pip1 <- susie_get_pip(fit1); pip2 <- susie_get_pip(fit2)
cs1  <- susie_get_cs(fit1, Xcorr = R); cs2  <- susie_get_cs(fit2, Xcorr = R)

cs_members <- function(cs, j) {
  if (is.null(cs$cs_index)) return(character(0))
  hit <- which(cs$cs_index == j)
  if (!length(hit)) return(character(0))
  snp_names[ cs$cs[[hit[1]]] ]
}

# --- 2) Top shared pairs (take top K) ---
K <- min(5L, nrow(pair))
top_pairs <- pair[order(-PP_shared)][1:K]

top_pairs[, `:=`(
  pip1_best = pip1[ match(best_snp, snp_names) ],
  pip2_best = pip2[ match(best_snp, snp_names) ],
  in_cs1    = mapply(function(j, s) s %in% cs_members(cs1, j), j, best_snp),
  in_cs2    = mapply(function(k, s) s %in% cs_members(cs2, k), k, best_snp),
  cs1_snp   = lapply(j, function(j) cs_members(cs1, j)),
  cs2_snp   = lapply(k, function(k) cs_members(cs2, k))
)]

# Jaccard of CS overlap (if both CS exist)
top_pairs[, jaccard := mapply(function(a,b){
  if (!length(a) || !length(b)) return(NA_real_)
  inter <- length(intersect(a,b)); uni <- length(unique(c(a,b)))
  if (uni == 0) NA_real_ else inter/uni
}, cs1_snp, cs2_snp)]

# Direction at best_snp (needs 'aln' from earlier alignment)
if (!"SNP" %in% names(aln)) aln[, SNP := snp_names]  # safeguard
dir_df <- aln[, .(SNP, beta_gwas, slope_aligned)]
top_pairs <- merge(top_pairs, dir_df, by.x = "best_snp", by.y = "SNP", all.x = TRUE)
top_pairs[, dir := ifelse(is.finite(beta_gwas) & is.finite(slope_aligned),
                          ifelse(sign(beta_gwas) == sign(slope_aligned), "concordant","discordant"),
                          NA_character_)]

# Keep just the tidy columns for the log
log_pairs <- top_pairs[, .(
  pair      = paste0("GWAS L", j, " ↔ sQTL L", k),
  PP_shared = signif(PP_shared, 3),
  best_snp,
  pip_gwas  = signif(pip1_best, 3),
  pip_sqtl  = signif(pip2_best, 3),
  in_cs1, in_cs2,
  jaccard   = ifelse(is.na(jaccard), NA, signif(jaccard, 3)),
  dir
)]

cat("\n# Top SuSiE×SuSiE shared component pairs\n")
print(log_pairs)

# --- 3) Per-trait component summary (lead SNP, PIP, CS size) ---
lead_of <- function(fit, pip) {
  L <- nrow(fit$alpha); out <- vector("list", L)
  for (j in seq_len(L)) {
    lead_idx <- which.max(fit$alpha[j, ])
    out[[j]] <- data.table(component = j,
                           lead_snp  = snp_names[lead_idx],
                           pip_lead  = pip[lead_idx])
  }
  rbindlist(out)
}
comp1 <- lead_of(fit1, pip1); comp2 <- lead_of(fit2, pip2)

cs_size <- function(cs, L) {
  out <- rep(NA_integer_, L)
  if (!is.null(cs$cs_index)) {
    for (i in seq_along(cs$cs_index)) out[cs$cs_index[i]] <- length(cs$cs[[i]])
  }
  out
}
comp1[, cs_size := cs_size(cs1, nrow(fit1$alpha))]
comp2[, cs_size := cs_size(cs2, nrow(fit2$alpha))]

cat("\n# Trait summaries\n")
cat("GWAS components:\n"); print(comp1)
cat("sQTL components:\n"); print(comp2)

# --- neat tables for the log (assumes: fit1, fit2, R, aln, z1, z2 exist) ---
library(data.table)
library(knitr)

# Names for SNP columns in alpha (prefer rsid, fallback to SNP, else index)
snp_names <- if ("rsid" %in% names(aln)) aln$rsid else if ("SNP" %in% names(aln)) aln$SNP else paste0("s", seq_len(ncol(R)))

# Posterior allocations per component (L x p) and PIPs
a1   <- fit1$alpha
a2   <- fit2$alpha
pip1 <- susie_get_pip(fit1)
pip2 <- susie_get_pip(fit2)

# Credible sets (index-based)
cs1 <- susie_get_cs(fit1, Xcorr = R)
cs2 <- susie_get_cs(fit2, Xcorr = R)

L1 <- nrow(a1); L2 <- nrow(a2)
stopifnot(L1 > 0, L2 > 0, length(snp_names) == ncol(a1), ncol(a1) == ncol(a2))

# Helper: size of CS j (or NA)
cs_size1 <- function(j) if (!is.null(cs1$cs) && length(cs1$cs) >= j && !is.null(cs1$cs[[j]])) length(cs1$cs[[j]]) else NA_integer_
cs_size2 <- function(k) if (!is.null(cs2$cs) && length(cs2$cs) >= k && !is.null(cs2$cs[[k]])) length(cs2$cs[[k]]) else NA_integer_

# Build pair table with PP_shared = sum_s alpha1[j,s]*alpha2[k,s]
pairs <- CJ(j = seq_len(L1), k = seq_len(L2))[, {
  w     <- a1[j, ] * a2[k, ]
  pp    <- sum(w, na.rm = TRUE)
  idx   <- which.max(w)
  snp   <- snp_names[idx]
  dir   <- if (!is.null(z1) && !is.null(z2) && length(z1) >= idx && length(z2) >= idx) {
    if (sign(z1[idx]) == sign(z2[idx])) "concordant" else "discordant"
  } else NA_character_
  # CS sizes and Jaccard
  csj <- if (!is.null(cs1$cs) && length(cs1$cs) >= j) cs1$cs[[j]] else NULL
  csk <- if (!is.null(cs2$cs) && length(cs2$cs) >= k) cs2$cs[[k]] else NULL
  jac <- if (!is.null(csj) && !is.null(csk) && length(csj) > 0 && length(csk) > 0) {
    inter <- length(intersect(csj, csk))
    union <- length(unique(c(csj, csk)))
    if (union > 0) inter/union else NA_real_
  } else NA_real_
  list(
    PP_shared = pp,
    best_snp  = snp,
    pip_gwas  = if (length(pip1) >= idx) pip1[idx] else NA_real_,
    pip_sqtl  = if (length(pip2) >= idx) pip2[idx] else NA_real_,
    cs1_size  = cs_size1(j),
    cs2_size  = cs_size2(k),
    jaccard   = jac,
    direction = dir
  )
}, by = .(j, k)]

setorder(pairs, -PP_shared)

# ---------- Table 1: Top shared component pairs ----------
cat("\n### Top SuSiE×SuSiE shared component pairs\n")
top_pairs <- pairs[1:min(15, .N)][
  , .(
    pair      = sprintf("GWAS L%d ↔ sQTL L%d", j, k),
    PP_shared = signif(PP_shared, 3),
    best_snp,
    pip_gwas  = signif(pip_gwas, 3),
    pip_sqtl  = signif(pip_sqtl, 3),
    cs1_size,
    cs2_size,
    jaccard   = ifelse(is.na(jaccard), NA, round(jaccard, 3)),
    direction
  )
]
print(kable(top_pairs, format = "markdown",
            caption = "Top SuSiE×SuSiE shared component pairs (PP_shared, best SNP, PIPs, CS sizes, Jaccard, direction)"))

# ---------- Table 2: GWAS components (SuSiE) ----------
cat("\n### GWAS components (SuSiE)\n")
gwas_comp <- rbindlist(lapply(seq_len(L1), function(j) {
  idx <- which.max(a1[j, ])
  data.table(
    component = j,
    lead_snp  = snp_names[idx],
    pip_lead  = signif(a1[j, idx], 3),
    cs_size   = cs_size1(j)
  )
}))
print(kable(gwas_comp, format = "markdown",
            caption = "GWAS components: lead SNP (by component posterior), lead PIP, and credible set size"))

# ---------- Table 3: sQTL components (SuSiE) ----------
cat("\n### sQTL components (SuSiE)\n")
sqtl_comp <- rbindlist(lapply(seq_len(L2), function(k) {
  idx <- which.max(a2[k, ])
  data.table(
    component = k,
    lead_snp  = snp_names[idx],
    pip_lead  = signif(a2[k, idx], 3),
    cs_size   = cs_size2(k)
  )
}))
print(kable(sqtl_comp, format = "markdown",
            caption = "sQTL components: lead SNP (by component posterior), lead PIP, and credible set size"))