SSDH_REPO <- getwd()   # repository root; run these scripts from there
setwd("/n/data1/bch/genetics/lee/shulin/heart_endo")

.libPaths("/n/data1/bch/genetics/lee/shulin/heart_endo/Rpackages_O2_upgraded/")
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(scales)
library(reshape2)
library(glmmTMB)
library(ggpubr)
library(rstatix)
library(Hmisc)
library(readxl)
library(mgcv)
library(ggeffects)

.libPaths("/home/shm283/shulin/softwares/anaconda3/envs/mutationalpatterns/lib/R/library")
library(MutationalPatterns)
ref_genome="BSgenome.Hsapiens.UCSC.hg19"
chr_orders=c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
library(ref_genome,character.only=TRUE)

source(file.path(SSDH_REPO, "analysis/utils.R"))

.libPaths("/n/data1/bch/genetics/lee/shulin/heart_endo/Rpackages_O2_upgraded/")
library(Hmisc)

# font
font = "Helvetica"
font_size = 7
axis_font_size = 6

# Custom function to format scientific notation with superscripts and handle zero
scientific_10 <- function(x) {
  sapply(x, function(v) {
    if (is.na(v)) {
      NA_character_
    } else if (v == 0) {
      "0"
    } else if (v < 0.05) {
      #parse(text = gsub("e", " %*% 10^", sprintf("%.2e", v)))
      lab <- sprintf("%.2e", v)
      parts <- strsplit(lab, "e")[[1]]
      
      mantissa <- parts[1]
      exponent <- as.integer(parts[2])
      
      parse(text = sprintf("'%s' %%*%% 10^%d", mantissa, exponent))
    } else {
      sprintf("%.2f", v)
    }
  })
}

scientific_10_plot <- function(x) {
  sapply(x, function(v) {
    if (is.na(v)) {
      NA_character_
    } else if (v == 0) {
      "0"
    } else if (v < 0.01) {
      parse(text = gsub("e", " %*% 10^", sprintf("%.2e", v)))
    } else {
      sprintf("%.2f", v)
    }
  })
}

condition_color_set <- c("Control" = "dodgerblue1", "AD" = "violetred3", "CTE" = "chocolate4")

# Meta data
metadata <- read.csv("./data/neuron_scan2/neuron_metadata_cell.csv")

# COSMIC signature 
snv_signatures <- readRDS("../ref/COSMIC_signatures/COSMIC_v3.2_SBS_GRCh37.rds")
indel_signatures <- readRDS("../ref/COSMIC_signatures/COSMIC_v3.2_ID_GRCh37.rds")

# Neuron signature (Michael A. Lodato et al., Science, 2018)
aging_signature_all <- readRDS("/n/data1/bch/genetics/lee/shulin/CTE/ref/aging_signature_all.rds")
aging_signature_all <- aging_signature_all[,1:3]
colnames(aging_signature_all) <- paste0("Signature ", c("A", "B", "C"))

genomic_factor <- 5.845

# ------------------------------------------------------------------------------
# Single-stranded damage burden (control, AD, CTE neurons)
single_strand_damage_burden <- readRDS("./results/ssDNA_mut_burden/neuron_single_strand_damage_burden.rds")
single_strand_damage_burden <- single_strand_damage_burden[single_strand_damage_burden$Diagnosis %in% c("Control", "AD", "CTE"),]

# Remove cells with too small callable regions (bottom 5%)
single_strand_damage_burden <- single_strand_damage_burden[single_strand_damage_burden$callable_region_size_filter,]

# Mixed effect model
ssSNV_lm <- glmmTMB(
  snv ~ Age + Diagnosis +
    offset(log(ssSNV_effective_callable_region_size)) +
    (1 | Case_ID),
  data = single_strand_damage_burden,
  family = nbinom2
)
summary(ssSNV_lm) -> ssSNV_lm_summary
p_AD_ssSNV <- ssSNV_lm_summary$coefficients$cond["DiagnosisAD", "Pr(>|z|)"]
p_CTE_ssSNV <- ssSNV_lm_summary$coefficients$cond["DiagnosisCTE", "Pr(>|z|)"]

pred <- ggpredict(ssSNV_lm, terms = "Age", 
                  condition = c(ssSNV_effective_callable_region_size = 1), 
                  type = "fixed")
rate_age_lm <- lm(predicted ~ x, data = pred)
ssSNV_age_slope <- coef(rate_age_lm)["x"]
ssSNV_age_intercept <- coef(rate_age_lm)["(Intercept)"]

ssIndel_lm <- glmmTMB(
  indel ~ Age + Diagnosis +
    offset(log(ssIndel_effective_callable_region_size)) +
    (1 | Case_ID),
  data = single_strand_damage_burden,
  family = nbinom2
)
summary(ssIndel_lm) -> ssIndel_lm_summary
p_AD_ssIndel <- ssIndel_lm_summary$coefficients$cond["DiagnosisAD", "Pr(>|z|)"]
p_CTE_ssIndel <- ssIndel_lm_summary$coefficients$cond["DiagnosisCTE", "Pr(>|z|)"]

pred <- ggpredict(ssIndel_lm, terms = "Age", 
                  condition = c(ssIndel_effective_callable_region_size = 1), 
                  type = "fixed")
rate_age_lm <- lm(predicted ~ x, data = pred)
ssIndel_age_slope <- coef(rate_age_lm)["x"]
ssIndel_age_intercept <- coef(rate_age_lm)["(Intercept)"]

# Plot
# SNV
ggplot(single_strand_damage_burden, 
       aes(x = Age, y = ssSNV_rate, color = Diagnosis)) +
  geom_point(size=0.5,shape=16) +
  geom_abline(intercept=ssSNV_age_intercept,slope=ssSNV_age_slope,
              size=0.5,color="lightskyblue") + 
  xlab("Age (yrs)") + ylab("ssDNA SNVs per bp") + 
  annotate(geom="text", x=0, y=2.5e-7*c(1, 0.925), 
           label = c(
             paste0("AD:~italic(P)==", scientific_10(p_AD_ssSNV)),
             paste0("CTE:~italic(P)==", scientific_10(p_CTE_ssSNV))
           ), 
           color = "black", family = font, size = 2, parse = T, hjust = 0) +
  scale_color_manual(values = condition_color_set) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(color = "black", size = font_size, family = font)
  ) +
  scale_y_continuous(labels = scientific_10_plot)
ggsave("./figures/ssDNA_mut_manuscript_figures/figure2/ssSNV_burden.pdf", width = 5, height = 4.5, units = "cm")

# Indel
ggplot(single_strand_damage_burden, 
       aes(x = Age, y = ssIndel_rate, color = Diagnosis)) +
  geom_point(size=0.5,shape=16) +
  geom_abline(intercept=ssIndel_age_intercept,slope=ssIndel_age_slope,
              size=0.5,color="lightskyblue") +
  xlab("Age (yrs)") + ylab("ssDNA Indels per bp") + 
  annotate(geom="text", x=0, y=2.5e-7*c(1, 0.925), 
           label = c(
             paste0("AD:~italic(P)==", scientific_10(p_AD_ssIndel)),
             paste0("CTE:~italic(P)==", scientific_10(p_CTE_ssIndel))
           ), 
           color = "black", family = font, size = 2, parse = T, hjust = 0) +
  scale_color_manual(values = condition_color_set) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(color = "black", size = font_size, family = font)
  ) +
  scale_y_continuous(labels = scientific_10_plot)
ggsave("./figures/ssDNA_mut_manuscript_figures/figure2/ssIndel_burden.pdf", width = 5, height = 4.5, units = "cm")

# ------------------------------------------------------------------------------
# Single-stranded damage spectrum (control, AD, CTE neurons)
single_strand_damage_mat <- readRDS("./results/ssDNA_mut_spectrum/neuron_single_strand_damage_mat_age_matched.rds")

# SNV
single_strand_SNV_mat_condition <- single_strand_damage_mat[["single_strand_SNV_mat_condition"]]
single_strand_SNV_mat_condition_diff <- single_strand_damage_mat[["single_strand_SNV_mat_condition_diff"]]
single_strand_SNV_mat_condition <- cbind(single_strand_SNV_mat_condition, single_strand_SNV_mat_condition_diff, aging_signature_all, 0)
colnames(single_strand_SNV_mat_condition)[dim(single_strand_SNV_mat_condition)[2]] <- "placeholder" # a placeholder for signature C (only for plot use)
colnames(single_strand_SNV_mat_condition)[1] <- "Control\n(Age-matched)" 
single_strand_SNV_mat_condition[single_strand_SNV_mat_condition < 0] <- 0

my.plot_96_profile_absolute(single_strand_SNV_mat_condition[, c("Control\n(Age-matched)", "AD", "CTE", "AD - Control", "CTE - Control", "placeholder")], 
                            condensed = T, facet_grid_scales = "free_x") + 
  theme_classic() +
  theme(
    strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    axis.text.x = element_text(color = "black", size = 2, family = font, angle = 90, vjust = 0.5),
    axis.ticks.x=element_blank(), 
    panel.spacing.x = unit(0.1, "lines"),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
  ) +
  labs(y = "Absolute contribution (ssDNA SNVs per bp)", x = "Context") +
  scale_y_continuous(labels = scientific_10_plot, n.breaks = 3)
ggsave("./figures/ssDNA_mut_manuscript_figures/figure2/ssSNV_spectrum.pdf", width = 9, height = 8, units = "cm")

my.plot_96_profile_absolute(single_strand_SNV_mat_condition[, c("Control\n(Age-matched)", "AD", "CTE", "AD - Control", "CTE - Control", "Signature C")],
                            condensed = T) + 
  theme_classic() +
  theme(
    strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    axis.text.x = element_text(color = "black", size = 2, family = font, angle = 90, vjust = 0.5),
    axis.ticks.x=element_blank(), 
    panel.spacing.x = unit(0.1, "lines"),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
  ) +
  labs(y = "Relative contribution", x = "Context") +
  scale_y_continuous(labels = scientific_10_plot)
ggsave("./figures/ssDNA_mut_manuscript_figures/figure2/SigC_spectrum.pdf", width = 9, height = 8, units = "cm")

# Indel
single_strand_Indel_mat_condition <- single_strand_damage_mat[["single_strand_Indel_mat_condition"]]
single_strand_Indel_mat_condition_diff <- single_strand_damage_mat[["single_strand_Indel_mat_condition_diff"]]
single_strand_Indel_mat_condition <- cbind(single_strand_Indel_mat_condition, single_strand_Indel_mat_condition_diff, 
                                           indel_signatures[, "ID4", drop = F], 0)
colnames(single_strand_Indel_mat_condition)[dim(single_strand_Indel_mat_condition)[2]] <- "placeholder" # a placeholder for ID4 (only for plot use)
colnames(single_strand_Indel_mat_condition)[1] <- "Control\n(Age-matched)" 
single_strand_Indel_mat_condition[single_strand_Indel_mat_condition < 0] <- 0

my.plot_indel_contexts(single_strand_Indel_mat_condition[, c("Control\n(Age-matched)", "AD", "CTE", "AD - Control", "CTE - Control", "placeholder")], 
                       condensed = T, same_y = T) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = 2, family = font, angle = 90, vjust = 0.5),
    axis.ticks.x = element_blank(), 
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    panel.spacing.x = unit(0.1, "lines"),
    legend.position = "none"
  ) +
  labs(y = "Absolute contribution (ssDNA Indels per bp)", x = "Context") +
  scale_y_continuous(labels = scientific_10_plot)
ggsave("./figures/ssDNA_mut_manuscript_figures/figure2/ssIndel_spectrum.pdf", width = 9, height = 8.38, units = "cm")

my.plot_indel_contexts(single_strand_Indel_mat_condition [, c("Control\n(Age-matched)", "AD", "CTE", "AD - Control", "CTE - Control", "ID4")], 
                       condensed = T, same_y = F) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = 2, family = font, angle = 90, vjust = 0.5),
    axis.ticks.x = element_blank(), 
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    panel.spacing.x = unit(0.1, "lines"),
    legend.position = "none"
  ) +
  labs(y = "Relative contribution", x = "Context") +
  scale_y_continuous(labels = scientific_10_plot)
ggsave("./figures/ssDNA_mut_manuscript_figures/figure2/ID4_spectrum.pdf", width = 9.09, height = 8.38, units = "cm")

# ------------------------------------------------------------------------------
# Single-stranded signature analysis (control, AD, CTE neurons)
# burden
single_strand_damage_burden <- readRDS("./results/ssDNA_mut_burden/neuron_single_strand_damage_burden.rds")
single_strand_damage_burden <- single_strand_damage_burden[single_strand_damage_burden$callable_region_size_filter,]
single_strand_damage_burden <- single_strand_damage_burden[single_strand_damage_burden$Diagnosis %in% c("Control", "AD", "CTE"),]

# spectrum by cell
single_strand_damage_mat_cell <- readRDS("./results/ssDNA_mut_spectrum/neuron_single_strand_damage_mat_cell.rds")
single_strand_SNV_mat_cell <- single_strand_damage_mat_cell$single_strand_SNV_mat_cell[, single_strand_damage_burden$Cell_ID]
single_strand_Indel_mat_cell <- single_strand_damage_mat_cell$single_strand_Indel_mat_cell[, single_strand_damage_burden$Cell_ID]

# Signature C
aging_signature_all <- aging_signature_all[row.names(single_strand_SNV_mat_cell),]

SNV_signature_fit <- fit_to_signatures(single_strand_SNV_mat_cell, aging_signature_all)
apply(SNV_signature_fit$contribution, 2, function(X){X/(sum(X)+1e-31)}) -> SNV_signature_fit
SNV_signature_fit <- data.frame(t(SNV_signature_fit))
SNV_signature_fit[,"Cell_ID"] <- row.names(SNV_signature_fit)
SNV_signature_fit <- right_join(single_strand_damage_burden, SNV_signature_fit)

SNV_signature_fit$Signature.C_tot_count <- SNV_signature_fit$Signature.C * SNV_signature_fit$snv
SNV_signature_fit$Signature.C_rate <- SNV_signature_fit$Signature.C_tot_count / SNV_signature_fit$ssSNV_effective_callable_region_size

# mixed effect model
ssSigC_lm <- glmmTMB(
  Signature.C_tot_count ~ Age + Diagnosis +
    offset(log(ssSNV_effective_callable_region_size)) +
    (1 | Case_ID),
  data = SNV_signature_fit,
  family = nbinom2
)
summary(ssSigC_lm) -> ssSigC_lm_summary
p_AD_ssSigC <- ssSigC_lm_summary$coefficients$cond["DiagnosisAD", "Pr(>|z|)"]
p_CTE_ssSigC <- ssSigC_lm_summary$coefficients$cond["DiagnosisCTE", "Pr(>|z|)"]

pred <- ggpredict(ssSigC_lm, terms = "Age", 
                  condition = c(ssSNV_effective_callable_region_size = 1), 
                  type = "fixed")
rate_age_lm <- lm(predicted ~ x, data = pred)
ssSigC_age_slope <- coef(rate_age_lm)["x"]
ssSigC_age_intercept <- coef(rate_age_lm)["(Intercept)"]

# plot
ggplot(SNV_signature_fit,aes(x=Age,y=Signature.C_rate,color=Diagnosis)) +
  geom_point(size=0.5,shape=16) +
  xlab("Age (yrs)") + ylab("ssDNA signature C contribution per bp") + 
  geom_abline(intercept=ssSigC_age_intercept,slope=ssSigC_age_slope,
              size=0.5,color="lightskyblue") + 
  annotate(geom="text", x=0, y=1.52e-7*c(1, 0.925), 
           label = c(
             paste0("AD:~italic(P)==", scientific_10(p_AD_ssSigC)),
             paste0("CTE:~italic(P)==", scientific_10(p_CTE_ssSigC))
           ), 
           color = "black", family = font, size = 2, parse = T, hjust = 0) +
  scale_color_manual(values = condition_color_set) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(color = "black", size = font_size, family = font)
  ) +
  scale_y_continuous(labels = scientific_10_plot)
ggsave("./figures/ssDNA_mut_manuscript_figures/figure2/ssSigC_burden.pdf", width = 5, height = 4.5, units = "cm")

# ID4
# Signature fitting
Indel_signature_fit <- fit_to_signatures(single_strand_Indel_mat_cell, indel_signatures)
apply(Indel_signature_fit$contribution, 2, function(X){X/(sum(X)+1e-31)}) -> Indel_signature_fit
Indel_signature_fit <- data.frame(t(Indel_signature_fit))
Indel_signature_fit[,"Cell_ID"] <- row.names(Indel_signature_fit)
Indel_signature_fit <- left_join(single_strand_damage_burden, Indel_signature_fit)

Indel_signature_fit$ID4_tot_count <- Indel_signature_fit$ID4 * Indel_signature_fit$indel
Indel_signature_fit$ID4_rate <- Indel_signature_fit$ID4_tot_count / Indel_signature_fit$ssIndel_effective_callable_region_size

# High ID4 cell
ID4_95_quantile <- quantile(Indel_signature_fit[Indel_signature_fit$Diagnosis == "Control", "ID4_rate"], 0.95)

# mixed effect model
ssID4_lm <- glmmTMB(
  ID4_tot_count ~ Age + Diagnosis +
    offset(log(ssIndel_effective_callable_region_size)) +
    (1 | Case_ID),
  data = Indel_signature_fit,
  family = nbinom2
)
summary(ssID4_lm) -> ssID4_lm_summary
p_AD_ssID4 <- ssID4_lm_summary$coefficients$cond["DiagnosisAD", "Pr(>|z|)"]
p_CTE_ssID4 <- ssID4_lm_summary$coefficients$cond["DiagnosisCTE", "Pr(>|z|)"]

pred <- ggpredict(ssID4_lm, terms = "Age", 
                  condition = c(ssIndel_effective_callable_region_size = 1), 
                  type = "fixed")
rate_age_lm <- lm(predicted ~ x, data = pred)
ssID4_age_slope <- coef(rate_age_lm)["x"]
ssID4_age_intercept <- coef(rate_age_lm)["(Intercept)"]

# plot
ggplot(Indel_signature_fit,aes(x=Age,y=ID4_rate,color=Diagnosis)) +
  geom_point(size=0.5,shape=16) +
  xlab("Age (yrs)") + ylab("ssDNA ID4 contribution per bp") + 
  geom_abline(intercept=ssID4_age_intercept,slope=ssID4_age_slope,
              size=0.5,color="lightskyblue") + 
  annotate(geom="text", x=0, y=2.1e-7*c(1, 0.925), 
           label = c(
             paste0("AD:~italic(P)==", scientific_10(p_AD_ssID4)),
             paste0("CTE:~italic(P)==", scientific_10(p_CTE_ssID4))
           ), 
           color = "black", family = font, size = 2, parse = T, hjust = 0) +
  scale_color_manual(values = condition_color_set) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(color = "black", size = font_size, family = font)
  ) +
  scale_y_continuous(labels = scientific_10_plot)
ggsave("./figures/ssDNA_mut_manuscript_figures/figure2/ssID4_burden.pdf", width = 5, height = 4.5, units = "cm")

# Correlation between signature C and ID4
inner_join(SNV_signature_fit, Indel_signature_fit) -> signature_fit
signature_fit$is_ID4_high <- ifelse(signature_fit$ID4_rate > ID4_95_quantile, "High", "Low")
signature_fit$is_ID4_high[signature_fit$Diagnosis == "Control"] <- "Control"

signature_fit <- signature_fit[signature_fit$Age > 35,] # Age-matched
signature_fit$is_ID4_high <- factor(signature_fit$is_ID4_high, 
                                    levels = c("Control", "Low", "High"),
                                    labels = c("Control cell\n(Age-matched)", "Low-ssDNA-ID4 cell\n(AD & CTE)", "High-ssDNA-ID4 cell\n(AD & CTE)"))

signature_fit$Diagnosis_plot <- factor(signature_fit$Diagnosis,
                                        levels = rev(c("Control", "AD", "CTE")),
                                        labels = rev(c("Control\n(Age-matched)", "AD", "CTE")))
ggplot(signature_fit,aes(x=Diagnosis_plot,y=ID4_rate,color=Diagnosis_plot)) +
  geom_jitter(size=0.5,shape=16, width = 0.3, height = 0) +
  labs(x = NULL, y = "ssDNA ID4 contribution per bp") + 
  geom_hline(yintercept = ID4_95_quantile, linewidth = 0.5, 
             color = "black", linetype = "dashed",) +
  scale_color_manual(values = rev(c("dodgerblue1",  "violetred3",  "chocolate4"))) +
  annotate(geom="text", x=3.3, y=ID4_95_quantile*1.5, 
           label = c("95th percentile of controls"), 
           family = font, size = 2, parse = F, hjust = 0) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(color = "black", size = font_size, family = font)
  ) +
  scale_y_continuous(labels = scientific_10_plot, limits = c(-3e-8, 2.2e-7)) +
  coord_flip()
ggsave("./figures/ssDNA_mut_manuscript_figures/figure2/ssID4_burden_group.pdf", width = 5, height = 4, units = "cm")

ssSigC_ssID4_lm <- lm(Signature.C_rate ~ ID4_rate, 
                      data = signature_fit[signature_fit$Diagnosis %in% c("AD", "CTE") & 
                                             signature_fit$is_ID4_high == "High-ID4 cell\n(AD & CTE)",])
ssSigC_ssID4_lm_summary <- summary(ssSigC_ssID4_lm)
ssSigC_ssID4_intercept <- ssSigC_ssID4_lm_summary$coefficients["(Intercept)", "Estimate"]
ssSigC_ssID4_slope <- ssSigC_ssID4_lm_summary$coefficients["ID4_rate", "Estimate"]
ssSigC_ssID4_pvalue <- ssSigC_ssID4_lm_summary$coefficients["ID4_rate", "Pr(>|t|)"]

ggplot(signature_fit[signature_fit$Diagnosis %in% c("AD", "CTE") & 
                       signature_fit$is_ID4_high == "High-ssDNA-ID4 cell\n(AD & CTE)",],
       aes(y = Signature.C_rate, x = ID4_rate)) +
  geom_point(size=0.5,shape=16, color = "darkred") +
  geom_abline(intercept=ssSigC_ssID4_intercept,
              slope=ssSigC_ssID4_slope,
              size=0.5,color="black") + 
  annotate(geom="text", x=0, y=1.4e-7,
           label = paste0("italic(P)==", scientific_10(ssSigC_ssID4_pvalue)), 
           color = "black", family = font, size = 2, parse = T, hjust = 0) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(color = "black", size = font_size, family = font)
  ) +
  scale_x_continuous(labels = scientific_10_plot) +
  scale_y_continuous(labels = scientific_10_plot) +
  labs(x = "ssDNA ID4 contribution per bp\n(High-ssDNA-ID4 cells only)", y = "ssDNA signature C contribution per bp")
ggsave("./figures/ssDNA_mut_manuscript_figures/figure2/ssSigC_ssID4_burden_correlation_high_ssID4.pdf", width = 4.5, height = 5, units = "cm")

wilcox.test(x = signature_fit[signature_fit$Diagnosis %in% c("AD", "CTE", "Control") & 
                           signature_fit$is_ID4_high == "High-ssDNA-ID4 cell\n(AD & CTE)", "Signature.C_rate"],
       y = signature_fit[signature_fit$Diagnosis %in% c("AD", "CTE", "Control") & 
                           signature_fit$is_ID4_high == "Low-ssDNA-ID4 cell\n(AD & CTE)", "Signature.C_rate"])
wilcox.test(x = signature_fit[signature_fit$Diagnosis %in% c("AD", "CTE", "Control") & 
                           signature_fit$is_ID4_high == "Control cell\n(Age-matched)", "Signature.C_rate"],
       y = signature_fit[signature_fit$Diagnosis %in% c("AD", "CTE", "Control") & 
                           signature_fit$is_ID4_high == "Low-ssDNA-ID4 cell\n(AD & CTE)", "Signature.C_rate"])

ggplot(signature_fit[signature_fit$Diagnosis %in% c("AD", "CTE", "Control"),],
       aes(y = Signature.C_rate, x = is_ID4_high, fill = is_ID4_high)) +
  geom_violin(alpha = 0.7, width = 0.5) +
  #geom_jitter(alpha = 0.5, size = 0.4) +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),   # mean ± 1 SD
    geom = "pointrange", size = 0.2
  ) +
  scale_fill_manual(values = c(as.character(condition_color_set["Control"]), "tomato", "darkred")) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("High-ssDNA-ID4 cell\n(AD & CTE)", "Low-ssDNA-ID4 cell\n(AD & CTE)"), 
                                        c("Low-ssDNA-ID4 cell\n(AD & CTE)", "Control cell\n(Age-matched)")), 
                     size = 2, step.increase = 0.1) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(color = "black", size = font_size, family = font)
  ) +
  scale_y_continuous(labels = scientific_10_plot) +
  labs(x = "         ssDNA ID4 contribution", y = "ssDNA signature C contribution per bp")
ggsave("./figures/ssDNA_mut_manuscript_figures/figure2/ssSigC_ssID4_burden_correlation.pdf", width = 4.5, height = 5, units = "cm")

# ------------------------------------------------------------------------------
# ssDNA damage and dsDNA mutations
# dsDNA spectrum
ds_indel_list <- read_excel("./data/neuron_scan2/neuron_ds_mutation_list.xlsx", sheet = "PTA Mutation Sites (Indel)")
ds_indel_list$End <- ds_indel_list$Position + nchar(ds_indel_list$`REF Base`) - 1
cell_id_list <- unique(ds_indel_list$`Cell ID`)
ds_indel_grange_list <- list()
for (cell_id in cell_id_list){
  ds_indel_list_cell <- ds_indel_list[ds_indel_list$`Cell ID` == cell_id,]
  ds_indel_grange <- GRanges(seqnames = paste0("chr", ds_indel_list_cell$Chromosome),
                             ranges = IRanges(start = ds_indel_list_cell$Position,
                                              end = ds_indel_list_cell$End),
                             ref = ds_indel_list_cell$`REF Base`,
                             alt = ds_indel_list_cell$`ALT Base`)
  GenomeInfoDb::genome(ds_indel_grange) = 'hg19'
  ds_indel_grange_list[[cell_id]] <- ds_indel_grange
}
ds_indel_grange_list <- GRangesList(ds_indel_grange_list)

# Signature analysis
indel <- get_indel_context(ds_indel_grange_list, ref_genome)
mut_mat_indel <- count_indel_contexts(indel)

dsIndel_signature_fit = fit_to_signatures(mut_mat_indel, indel_signatures)
apply(dsIndel_signature_fit$contribution, 2, function(X){X/(sum(X)+1e-31)}) -> dsIndel_signature_fit
dsIndel_signature_fit <- data.frame(t(dsIndel_signature_fit))
dsIndel_signature_fit[,"Cell_ID"] <- rownames(dsIndel_signature_fit)

dsIndel_signature_fit <- left_join(metadata, dsIndel_signature_fit, by = "Cell_ID")
dsIndel_signature_fit$ds_ID4_rate <- dsIndel_signature_fit$ID4 * dsIndel_signature_fit$dsIndel_rate

# ssDNA signature analysis
Indel_signature_fit <- fit_to_signatures(single_strand_Indel_mat_cell, indel_signatures)
apply(Indel_signature_fit$contribution, 2, function(X){X/(sum(X)+1e-31)}) -> Indel_signature_fit
Indel_signature_fit <- data.frame(t(Indel_signature_fit))
Indel_signature_fit[,"Cell_ID"] <- row.names(Indel_signature_fit)
Indel_signature_fit <- left_join(single_strand_damage_burden, Indel_signature_fit)

Indel_signature_fit$ID4_tot_count <- Indel_signature_fit$ID4 * Indel_signature_fit$indel
Indel_signature_fit$ID4_rate <- Indel_signature_fit$ID4_tot_count / Indel_signature_fit$ssIndel_effective_callable_region_size

# combine ssDNA and dsDNA signature results
ds_ss_Indel_signature_fit <- inner_join(dsIndel_signature_fit, Indel_signature_fit, by = "Cell_ID", suffix = c("_ds", "_ss"))
ds_ss_Indel_signature_fit <- ds_ss_Indel_signature_fit[!is.na(ds_ss_Indel_signature_fit$ID4_ds),]

ds_ss_Indel_signature_fit <- ds_ss_Indel_signature_fit[ds_ss_Indel_signature_fit$Diagnosis_ds %in% c("AD", "CTE"),]
ds_ss_Indel_signature_fit$ssID4_rate_group <- ntile(ds_ss_Indel_signature_fit$ID4_rate, n = 10)
ds_ss_Indel_signature_fit_summary <- 
  ds_ss_Indel_signature_fit %>% group_by(ssID4_rate_group, Diagnosis_ds) %>% 
  summarise("mean_dsID4" = mean(ds_ID4_rate),"mean_ssID4" = mean(ID4_rate))

ggplot(ds_ss_Indel_signature_fit_summary, 
       aes(ssID4_rate_group, mean_dsID4/genomic_factor/1e9, color = Diagnosis_ds)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = condition_color_set) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(color = "black", size = font_size, family = font)
  ) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(labels = scientific_10_plot) +
  labs(x = "ssDNA ID4 contribution group", y = "dsDNA ID4 contribution per bp")
ggsave("./figures/ssDNA_mut_manuscript_figures/figure2/ssID4_dsID4_correlation.pdf", width = 6, height = 5, units = "cm")

