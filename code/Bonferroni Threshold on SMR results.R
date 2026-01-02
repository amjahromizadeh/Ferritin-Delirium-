library(data.table)

# ===========================
# INPUTS
# ===========================
INPUT_FILES <- c("Delirium and Blood.zip")

# Optional: write outputs (leave NULL if you only want objects in memory)
OUT_DIR <- NULL
# OUT_DIR <- "smr_bonf_outputs"

# Canonical filters (match colleague)
P_HEIDI_MIN <- 0.01
NSNP_MIN    <- 3

# ===========================
# Helpers 
# ===========================
infer_trait_omics <- function(filename) {
  trait <- NA_character_
  if (grepl("^Ferritin", filename, ignore.case = FALSE)) trait <- "Ferritin"
  if (grepl("^Delirium", filename, ignore.case = FALSE)) trait <- "Delirium"
  
  omics <- regmatches(filename, regexpr("(eQTL|mQTL|sQTL|pQTL)", filename, perl = TRUE))
  if (length(omics) == 0) omics <- NA_character_
  
  list(trait = trait, omics = omics)
}

assign_tissue_from_qtl <- function(qtl_name) {
  s <- tolower(trimws(as.character(qtl_name)))
  if (is.na(s) || s == "") return(NA_character_)
  if (grepl("pituitary", s)) return("Pituitary")
  if (grepl("kidney", s) && grepl("cortex", s)) return("KidneyCortex")
  if (grepl("spleen", s)) return("Spleen")
  if (grepl("liver", s)) return("Liver")
  if (grepl("brain", s) || grepl("brainmeta", s)) return("Brain")
  if (grepl("blood", s) || grepl("whole_blood", s) || grepl("eqtlgen", s) ||
      grepl("mcrae", s) || grepl("^pqtl_", s)) return("Blood")
  NA_character_
}

standardize_cols_minimal <- function(dt) {
  # Probe ID
  if ("probeID" %in% names(dt) && !"probe_id" %in% names(dt)) setnames(dt, "probeID", "probe_id")
  if ("probe_id" %in% names(dt) && !"probe_id" %in% names(dt)) setnames(dt, "probe_id", "probe_id") # noop
  
  # P-values
  if ("p_SMR" %in% names(dt) && !"p_smr" %in% names(dt)) setnames(dt, "p_SMR", "p_smr")
  if ("p_HEIDI" %in% names(dt) && !"p_heidi" %in% names(dt)) setnames(dt, "p_HEIDI", "p_heidi")
  
  # nsnp HEIDI (many aliases)
  nsnp_old <- intersect(names(dt), c("nsnp_HEIDI", "Nsnp(HEIDI)", "Nsnp", "nsnp", "Nsnp(heidi)"))
  if (length(nsnp_old) > 0 && !"nsnp_heidi" %in% names(dt)) {
    setnames(dt, nsnp_old[1], "nsnp_heidi")
  } else if (!"nsnp_heidi" %in% names(dt)) {
    dt[, nsnp_heidi := NA_real_]
  }
  
  # qtl_name (must exist for panel-wise BF)
  if (!"qtl_name" %in% names(dt)) dt[, qtl_name := NA_character_]
  
  # gene/top SNP naming (best-effort)
  if ("Gene" %in% names(dt) && !"gene" %in% names(dt)) setnames(dt, "Gene", "gene")
  if ("topSNP" %in% names(dt) && !"top_snp" %in% names(dt)) setnames(dt, "topSNP", "top_snp")
  
  # coerce numeric
  for (cc in c("p_smr", "p_heidi", "nsnp_heidi")) {
    if (cc %in% names(dt)) dt[[cc]] <- suppressWarnings(as.numeric(dt[[cc]]))
  }
  
  dt
}

read_smr_file <- function(path) {
  fn <- basename(path)
  meta <- infer_trait_omics(fn)
  
  if (grepl("\\.zip$", path, ignore.case = TRUE)) {
    tmp <- tempfile("smr_zip_")
    dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
    unzip(path, exdir = tmp)
    
    inner <- list.files(tmp, pattern = "\\.(tsv|txt|csv)$", full.names = TRUE, ignore.case = TRUE)
    if (length(inner) == 0) stop(sprintf("ZIP has no TSV/TXT/CSV: %s", fn))
    
    dt_list <- lapply(inner, fread)
    dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  } else {
    dt <- fread(path)
  }
  
  dt <- standardize_cols_minimal(dt)
  dt[, source_file := fn]
  
  # Attach trait/omics if inferable; otherwise leave for user to fill manually
  dt[, trait := meta$trait]
  dt[, omics := meta$omics]
  
  # Tissue derived from qtl_name 
  dt[, tissue := vapply(qtl_name, assign_tissue_from_qtl, character(1))]
  
  dt
}

# ===========================
# Load and combine
# ===========================
dt_all <- rbindlist(lapply(INPUT_FILES, read_smr_file), use.names = TRUE, fill = TRUE)

# Guardrails: these are required for panel-wise BF
req <- c("probe_id", "p_smr", "p_heidi", "nsnp_heidi", "qtl_name", "trait", "omics", "tissue")
missing_req <- setdiff(req, names(dt_all))
if (length(missing_req) > 0) stop(sprintf("Missing required columns: %s", paste(missing_req, collapse = ", ")))

# ===========================
# (1) Overall Bonferroni (ALL PANELS combined) 
# ===========================
M_overall   <- uniqueN(dt_all$probe_id)
alpha_all   <- 0.05 / M_overall

overall_summary <- data.table(
  scope         = "all panels combined",
  probes_tested = M_overall,
  alpha_bonf    = alpha_all
)

overall_summary[, n_pass_bonf_heidi := dt_all[
  p_smr <= alpha_all & p_heidi >= P_HEIDI_MIN & nsnp_heidi >= NSNP_MIN, .N
]]

# ===========================
# (2) Panel-wise Bonferroni per trait × tissue × omics × qtl_name 
# ===========================
dt_all[, M_tested_unique_probes := uniqueN(probe_id), by = .(trait, tissue, omics, qtl_name)]
dt_all[, p_bonf_panel := ifelse(M_tested_unique_probes > 0, 0.05 / M_tested_unique_probes, NA_real_)]

dt_all[, pass_bonf_heidi := !is.na(p_smr) & !is.na(p_heidi) & !is.na(nsnp_heidi) &
         (p_smr <= p_bonf_panel) & (p_heidi >= P_HEIDI_MIN) & (nsnp_heidi >= NSNP_MIN)]

passed_dt <- dt_all[pass_bonf_heidi == TRUE]

# Panel summary
panel_summary <- dt_all[, .(
  M_tested_unique_probes = uniqueN(probe_id),
  p_bonf_panel           = ifelse(uniqueN(probe_id) > 0, 0.05 / uniqueN(probe_id), NA_real_),
  n_pass_probes          = uniqueN(probe_id[pass_bonf_heidi == TRUE]),
  n_pass_unique_genes    = if ("gene" %in% names(dt_all)) uniqueN(gene[pass_bonf_heidi == TRUE]) else NA_integer_,
  n_pass_unique_probe_id = uniqueN(probe_id[pass_bonf_heidi == TRUE]),
  min_p_smr_pass         = if (any(pass_bonf_heidi == TRUE, na.rm = TRUE)) min(p_smr[pass_bonf_heidi == TRUE], na.rm = TRUE) else NA_real_,
  top_gene_pass          = if ("gene" %in% names(dt_all) && any(pass_bonf_heidi == TRUE, na.rm = TRUE)) gene[pass_bonf_heidi == TRUE][which.min(p_smr[pass_bonf_heidi == TRUE])] else NA_character_,
  top_probe_pass         = if (any(pass_bonf_heidi == TRUE, na.rm = TRUE)) probe_id[pass_bonf_heidi == TRUE][which.min(p_smr[pass_bonf_heidi == TRUE])] else NA_character_,
  top_snp_pass           = if ("top_snp" %in% names(dt_all) && any(pass_bonf_heidi == TRUE, na.rm = TRUE)) top_snp[pass_bonf_heidi == TRUE][which.min(p_smr[pass_bonf_heidi == TRUE])] else NA_character_
), by = .(trait, tissue, omics, qtl_name)][order(trait, tissue, omics, qtl_name)]

# ===========================
# (3) Extra “strict” table, but made consistent:
#     p_smr <= 1e-8 AND HEIDI-pass AND nsnp>=3 (independent of Bonferroni)
# ===========================
sig_strict <- dt_all[p_smr <= 1e-8 & p_heidi >= P_HEIDI_MIN & nsnp_heidi >= NSNP_MIN]
sig_strict_summary <- sig_strict[, .(n_hits = .N), by = .(trait, omics, tissue, qtl_name)][order(-n_hits)]

# ===========================
# (4) Collapsed summary across tissue × omics (matches colleague’s collapsed outputs)
# ===========================
collapsed_counts <- panel_summary[, .(
  M_tested_unique_probes = sum(M_tested_unique_probes, na.rm = TRUE),
  n_pass_probes          = sum(n_pass_probes, na.rm = TRUE)
), by = .(trait, tissue, omics)]

collapsed_gene_probe <- passed_dt[, .(
  n_pass_unique_genes    = if ("gene" %in% names(passed_dt)) uniqueN(gene) else NA_integer_,
  n_pass_unique_probe_id = uniqueN(probe_id)
), by = .(trait, tissue, omics)]

collapsed_summary <- merge(collapsed_counts, collapsed_gene_probe,
                           by = c("trait", "tissue", "omics"), all.x = TRUE)[order(trait, tissue, omics)]

# ===========================
# write outputs 
# ===========================
if (!is.null(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  fwrite(overall_summary,   file.path(OUT_DIR, "BF_overall_summary.tsv"), sep = "\t")
  fwrite(panel_summary,     file.path(OUT_DIR, "BF_summary_by_panel.tsv"), sep = "\t")
  fwrite(collapsed_summary, file.path(OUT_DIR, "BF_counts_trait_tissue_omics_collapsed.tsv"), sep = "\t")
  fwrite(passed_dt,         file.path(OUT_DIR, "BF_HEIDI_passed_hits.tsv"), sep = "\t")
  fwrite(sig_strict,        file.path(OUT_DIR, "HEIDI_passed_pSMR_1e-8_hits.tsv"), sep = "\t")
  fwrite(sig_strict_summary,file.path(OUT_DIR, "HEIDI_passed_pSMR_1e-8_counts_by_panel.tsv"), sep = "\t")
}


