library(data.table)
library(TwoSampleMR)
library(MRPRESSO)

exposure_raw <- fread("Ferritin_AF0p005.mr_ready.tsv.gz")
outcome_raw  <- fread("Delirium_AF0p005.mr_ready.tsv.gz") 

# Map to TwoSampleMR expected names
exp_dat <- copy(exposure_raw)
setnames(exp_dat, c("effect_allele","other_allele","beta","se","eaf","pval","samplesize"),
         c("effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure",
           "eaf.exposure","pval.exposure","samplesize.exposure"))
exp_dat[, exposure := "Ferritin"]

# Outcome mapping
out_dat <- copy(outcome_raw)
setnames(out_dat, c("effect_allele","other_allele","beta","se","eaf","pval","samplesize"),
         c("effect_allele.outcome","other_allele.outcome","beta.outcome","se.outcome",
           "eaf.outcome","pval.outcome","samplesize.outcome"))
out_dat[, outcome := "Delirium"]

# Start at p<5e-8; if <3 SNPs remain, relax to p<5e-6 as pre-specified sensitivity.
exp_gws <- exp_dat[pval.exposure < 5e-8]
if (nrow(exp_gws) < 3) exp_gws <- exp_dat[pval.exposure < 5e-6]

#Here, (OPENGWAS_JWT = "") must include a token from https://api.opengwas.io
Sys.setenv(OPENGWAS_JWT = "")

exp_gws <- clump_data(
  exp_gws,
  clump_kb = 10000,
  clump_r2 = 0.001,
  pop = "EUR")

# Harmonisation
exp_fmt <- format_data(
  as.data.frame(exp_gws), type = "exposure", snp_col = "SNP",
  beta_col = "beta.exposure", se_col = "se.exposure",
  eaf_col  = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col  = "other_allele.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  phenotype_col = "exposure"
)

# Keep only outcome rows at IV SNPs, then format
outcome_iv <- out_dat[SNP %in% exp_gws$SNP]
out_fmt <- format_data(
  as.data.frame(outcome_iv), type = "outcome", snp_col = "SNP",
  beta_col = "beta.outcome", se_col = "se.outcome",
  eaf_col  = "eaf.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col  = "other_allele.outcome",
  pval_col = "pval.outcome",
  samplesize_col = "samplesize.outcome",
  phenotype_col = "outcome"
)

exp_fmt$exposure <- "Ferritin"
out_fmt$outcome  <- "Delirium"

dat_h <- harmonise_data(exp_fmt, out_fmt, action = 2) 

# F-stat per SNP ~ (beta.exposure^2 / se.exposure^2); mean F as a quick check
dat_h$F_exposure <- (dat_h$beta.exposure^2) / (dat_h$se.exposure^2)
mean_F <- mean(dat_h$F_exposure, na.rm=TRUE); mean_F

# PRIMARY MR + SENSITIVITIES 
mr_main <- mr(dat_h, method_list = c("mr_ivw", "mr_ivw_mre",
                                     "mr_egger_regression", "mr_weighted_median"))
het       <- mr_heterogeneity(dat_h)               # Cochran's Q
pleio     <- mr_pleiotropy_test(dat_h)             # Egger intercept
loo <- mr_leaveoneout(dat_h)

# Plot 
p_loo <- mr_leaveoneout_plot(loo)
print(p_loo[[1]])
steiger   <- directionality_test(dat_h)            # Steiger

# MR-PRESSO (global test + outliers + distortion)
set.seed(1)
mrp <- mr_presso(BetaOutcome = "beta.outcome",
                 BetaExposure = "beta.exposure",
                 SdOutcome = "se.outcome",
                 SdExposure = "se.exposure",
                 OUTLIERtest = TRUE,
                 DISTORTIONtest = TRUE,
                 data = as.data.frame(dat_h),
                 NbDistribution = 1000,  SignifThreshold = 0.05)

# OUTPUTS
print(mr_main)
print(het)
print(pleio)
print(steiger)
print(mrp)

library(knitr)

# print markdown table
kable(mr_main, format = "markdown") 

#-------------------------------------------------------------------------------
library(ggplot2)

# 1) Pull rows from TwoSampleMR results
get_row <- function(d, pattern) d[grepl(pattern, d$method), ][1, ]

ivw_re <- get_row(mr_main, "Inverse variance weighted \\(multiplicative random effects\\)|Inverse variance weighted \\(random effects\\)")
egger  <- get_row(mr_main, "^MR Egger")
wm     <- get_row(mr_main, "^Weighted median")

# 2) Pull MR-PRESSO outlier-corrected estimate (if available)
presso_df <- tryCatch(as.data.table(mrp$`Main MR results`), error = function(e) NULL)
oc <- if (!is.null(presso_df)) presso_df[grepl("Outlier", `MR Analysis`, ignore.case = TRUE)] else NULL

# 3) Build a tidy table of estimates on the log-OR scale
est <- rbindlist(list(
  data.table(Method = "IVW (RE)",                               beta = ivw_re$b,          se = ivw_re$se),
  if (!is.null(oc) && nrow(oc)) data.table(Method = "IVW (MR-PRESSO outlier-corrected)",
                                           beta = as.numeric(oc$`Causal Estimate`),
                                           se   = as.numeric(oc$Sd)) else NULL,
  data.table(Method = "MR-Egger",                               beta = egger$b,           se = egger$se),
  data.table(Method = "Weighted median",                        beta = wm$b,              se = wm$se)
), use.names = TRUE, fill = TRUE)

# 4) Convert to ORs and 95% CIs
est[, `:=`(
  OR  = exp(beta),
  LCI = exp(beta - 1.96 * se),
  UCI = exp(beta + 1.96 * se)
)]
est[, Method := factor(Method, levels = c("IVW (RE)", "IVW (MR-PRESSO outlier-corrected)", "MR-Egger", "Weighted median"))]

x_min <- floor(min(est$LCI, na.rm = TRUE) * 100) / 100
x_max <- ceiling(max(est$UCI, na.rm = TRUE) * 100) / 100

p_mr1 <- ggplot(est, aes(x = OR, y = Method)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI), height = 0.15) +
  scale_x_continuous(limits = c(x_min, x_max),
                     expand = expansion(mult = c(0.02, 0.08))) +
  coord_cartesian(clip = "off") +                    # <- don't clip at panel edge
  labs(x = "Odds ratio (per SD higher Ferritin)", y = NULL,
       title = "Figure MR1. Causal estimates across MR estimators") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 8, r = 24, b = 8, l = 8)  # extra right margin
  )

print(p_mr1)
ggsave("Figure_MR1_forest.png", p_mr1,
       width = 7.5, height = 4.0, dpi = 300, limitsize = FALSE)
ggsave("Figure_MR1_forest.pdf",  p_mr1,
       width = 7.5, height = 4.0, useDingbats = FALSE)
print(p_mr1)
ggsave("Figure_MR1_forest.png", p_mr1,
       width = 7.5, height = 4.0, dpi = 300, limitsize = FALSE)
ggsave("Figure_MR1_forest.pdf",  p_mr1, width = 6.0, height = 3.5)
#-------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(scales)
})

# Pull MR-PRESSO outlier rsIDs 
get_presso_outliers <- function(mrp_obj){
  tryCatch({
    ot <- mrp_obj[["MR-PRESSO results"]][["Outlier Test"]][["Outliers"]]
    if (is.null(ot) || !is.data.frame(ot) || nrow(ot) == 0) return(character(0))
    pick <- intersect(c("SNP","rsid","Name","Outlier"), names(ot))
    if (length(pick) == 0) return(character(0))
    unique(as.character(ot[[pick[1]]]))
  }, error = function(e) character(0))
}

# pick one slope per method (prefer IVW random-effects if present)
pick_b <- function(df, exact_names, fallback_pattern = NULL){
  hits <- which(df$method %in% exact_names)
  if (length(hits) > 0) return(df$b[hits[1]])
  if (!is.null(fallback_pattern)) {
    hits <- which(grepl(fallback_pattern, df$method))
    if (length(hits) > 0) return(df$b[hits[1]])
  }
  NA_real_
}

b_ivw <- pick_b(
  mr_main,
  exact_names = c("Inverse variance weighted (multiplicative random effects)",
                  "Inverse variance weighted (random effects)"),
  fallback_pattern = "^Inverse variance weighted$"
)

b_egger <- pick_b(mr_main, exact_names = c("MR Egger"), fallback_pattern = "Egger")
b_wmed  <- pick_b(mr_main, exact_names = c("Weighted median"), fallback_pattern = "Weighted median")

egger_int <- if (!is.null(pleio$egger_intercept)) pleio$egger_intercept[1] else NA_real_
if (is.na(egger_int)) egger_int <- 0  # safe default if Egger intercept not returned

# one row per line to draw
line_df <- data.frame(
  method    = c("IVW (RE)", "MR-Egger", "Weighted median"),
  slope     = c(b_ivw,      b_egger,    b_wmed),
  intercept = c(0,          egger_int,  0),
  linetype  = c("solid",    "dashed",   "dotdash"),
  stringsAsFactors = FALSE
)

presso_outliers <- get_presso_outliers(mrp)

# ======================================================================
# S-MR1. MR scatter (with IVW / MR-Egger / weighted-median lines)
# ======================================================================
scatter_df <- as.data.frame(dat_h)
scatter_df$group <- ifelse(scatter_df$SNP %in% presso_outliers,
                           "PRESSO outlier", "Instrument")

line_df <- tibble(
  method    = c("IVW (RE)", "MR-Egger", "Weighted median"),
  slope     = c(b_ivw,      b_egger,    b_wmed),
  intercept = c(0,          egger_int,  0),
  linetype  = c("solid",    "dashed",   "dotdash")
)

p_scatter <- ggplot(scatter_df, aes(beta.exposure, beta.outcome)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(aes(shape = group), size = 2.4, alpha = 0.9) +
  geom_abline(data = line_df,
              aes(slope = slope, intercept = intercept, linetype = method),
              linewidth = 0.9, show.legend = TRUE) +
  scale_shape_manual(values = c("Instrument" = 16, "PRESSO outlier" = 17)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
  labs(title = "Supplementary Figure S-MR1. MR scatter plot",
       x = "SNP effect on Ferritin (beta)",
       y = "SNP effect on Delirium (beta)",
       linetype = "Estimator", shape = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(10, 70, 10, 10)) +
  coord_cartesian(clip = "off")

if (length(presso_outliers) > 0) {
  p_scatter <- p_scatter +
    ggrepel::geom_text_repel(
      data = subset(scatter_df, SNP %in% presso_outliers),
      aes(label = SNP), size = 3, max.overlaps = Inf
    )
}

ggsave("Supp_Figure_S-MR1_scatter.png", p_scatter,
       width = 7.5, height = 5.0, dpi = 300, bg = "white", limitsize = FALSE)

# ======================================================================
# S-MR2. Funnel plot (per-SNP ratio estimates)
# ======================================================================
single <- mr_singlesnp(dat_h)
p_funnel <- mr_funnel_plot(single)[[1]] +
  labs(title = "Supplementary Figure S-MR2. Funnel plot",
       x = "Causal estimate (Wald ratio per SNP)",
       y = "SE of ratio") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave("Supp_Figure_S-MR2_funnel.png", p_funnel,
       width = 7.0, height = 5.0, dpi = 300, bg = "white", limitsize = FALSE)

# ======================================================================
# S-MR3. Leave-one-out influence plot
# ======================================================================
loo <- mr_leaveoneout(dat_h)
p_loo <- mr_leaveoneout_plot(loo)[[1]] +
  labs(title = "Supplementary Figure S-MR3. Leave-one-out influence plot",
       y = "IVW estimate") +
  theme_minimal(base_size = 12)

ggsave("Supp_Figure_S-MR3_leaveoneout.png", p_loo,
       width = 7.0, height = 5.0, dpi = 300, bg = "white", limitsize = FALSE)

# ======================================================================
# S-MR4. MR-PRESSO outlier map (studentized residuals vs leverage)
# ======================================================================
ivw_fit <- lm(beta.outcome ~ beta.exposure,
              weights = 1/(se.outcome^2),
              data = scatter_df)

press_df <- scatter_df %>%
  mutate(leverage = hatvalues(ivw_fit),
         stud_res = rstudent(ivw_fit),
         outlier  = SNP %in% presso_outliers)

# Pull MR-PRESSO p-values (best-effort; shows NA if structure differs)
safe_get <- function(x, path, default = NA){
  tryCatch({ for (nm in path) x <- x[[nm]]; x }, error = function(e) default)
}
p_glob <- safe_get(mrp, c("MR-PRESSO results","Global Test","Pvalue"))
p_dist <- safe_get(mrp, c("MR-PRESSO results","Distortion Test","Pvalue"))

subtxt <- sprintf("Global p = %s; Distortion p = %s; Outliers = %d",
                  ifelse(is.na(p_glob), "NA", formatC(p_glob, format = "e", digits = 2)),
                  ifelse(is.na(p_dist), "NA", formatC(p_dist, format = "e", digits = 2)),
                  length(presso_outliers))

p_presso <- ggplot(press_df, aes(leverage, stud_res, colour = outlier)) +
  geom_hline(yintercept = c(-3, 0, 3),
             linetype = c("dotted", "solid", "dotted"), colour = "grey70") +
  geom_point(size = 2.3, alpha = 0.9) +
  scale_color_manual(values = c("FALSE" = "grey30", "TRUE" = "firebrick"),
                     labels = c("Instrument", "PRESSO outlier")) +
  labs(title = "Supplementary Figure S-MR4. MR-PRESSO outlier map",
       subtitle = subtxt,
       x = "Leverage (hat values from IVW fit)",
       y = "Studentized residuals",
       colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.margin = margin(10, 70, 10, 10)) +
  coord_cartesian(clip = "off")

if (length(presso_outliers) > 0) {
  p_presso <- p_presso +
    ggrepel::geom_text_repel(
      data = subset(press_df, outlier),
      aes(label = SNP), size = 3, max.overlaps = Inf
    )
}

ggsave("Supp_Figure_S-MR4_presso_map.png", p_presso,
       width = 7.2, height = 5.2, dpi = 300, bg = "white", limitsize = FALSE)
