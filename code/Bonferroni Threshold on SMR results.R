library(data.table)

# Read SMR portal results
dat <- fread("Delirium and Blood.zip")


n_probes <- length(unique(dat$probeID))
alpha_bonf <- 0.05 / n_probes 

#The subset of SMR results that survive Bonferroni (per tissue) and pass HEIDI
sig_bonf <- dat[p_SMR <= alpha_bonf & p_HEIDI >= 0.01]

#-------------------------------------------------------------------------------
# To get only p heidi >= 0.01 and p smr <= 1.e-8 regardless of Bonferroni
sig_strict <- sig_bonf[p_SMR <= 1e-8 & p_HEIDI >= 0.01]
sig_strict_summary <- sig_strict[,.(n_hits = .N),by = qtl_name][order(-n_hits)]

#-------------------------------------------------------------------------------
dt <- copy(dat)
setDT(dt)
# normalize HEIDI SNP-count column name
setnames(dt, old = intersect(names(dt), c("nsnp_HEIDI","Nsnp(HEIDI)","Nsnp")),
         new = "Nsnp", skip_absent = TRUE)

#Overall (all panels)
M_overall    <- uniqueN(dt$probeID)
alpha_bonf   <- 0.05 / M_overall
n_bonf_pass  <- dt[ p_SMR <= alpha_bonf & p_HEIDI >= 0.01 & Nsnp >= 3, .N ]

# how many probes were tested in total, 
# the corresponding Bonferroni threshold (0.05 / total probes), 
# and how many hits pass SMR+HEIDI under that threshold.
overall_summary <- data.table(
  scope          = "all panels combined",
  probes_tested  = M_overall,
  alpha_bonf     = alpha_bonf,
  n_bonf_pass    = n_bonf_pass
)
overall_summary

# Per tissue: probes tested (Bonferroni denominator),
# the per-tissue Bonferroni alpha,
# and the count of SMR+HEIDI Bonferroni-significant hits in that tissue.
per_tissue <- dt[, {
  M  <- uniqueN(probeID)
  a  <- 0.05 / M
  n  <- sum(p_SMR <= a & p_HEIDI >= 0.01 & Nsnp >= 3)
  .(probes_tested = M, alpha_bonf = a, n_bonf_pass = n)
}, by = qtl_name][order(qtl_name)]

#-------------------------------------------------------------------------------


