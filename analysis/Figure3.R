SSDH_REPO <- getwd()   # repository root; run these scripts from there
setwd("/n/data1/bch/genetics/lee/shulin/heart_endo")

.libPaths("/n/data1/bch/genetics/lee/shulin/heart_endo/Rpackages_O2_upgraded/")
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(scales)
library(reshape2)
library(glmmTMB)
library(ggpubr)
library(ggbreak)
library(ggsci)
library(ggeffects)
library(lme4)
library(lmerTest)

.libPaths("/home/shm283/shulin/softwares/anaconda3/envs/mutationalpatterns/lib/R/library")
library(MutationalPatterns)
ref_genome="BSgenome.Hsapiens.UCSC.hg19"
chr_orders=c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
library(ref_genome,character.only=TRUE)

source(file.path(SSDH_REPO, "analysis/utils.R"))

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
      parse(text = gsub("e", " %*% 10^", sprintf("%.2e", v)))
    } else {
      sprintf("%.2f", v)
    }
  })
}

condition_color_set <- c("Non-smoker w/o ASCVD" = "lightskyblue", "Smoker w/o ASCVD" = "dodgerblue4",
                         "Smoker w/ ASCVD" = "brown4")

genomic_factor <- 5.845

# ------------------------------------------------------------------------------
# Single-stranded damage burden (CVEC)
single_strand_damage_burden <- readRDS("./results/ssDNA_mut_burden/EC_single_strand_damage_burden.rds")
single_strand_damage_burden$Diagnosis <- str_c(single_strand_damage_burden$Smoking, "_", single_strand_damage_burden$ASCVD)

single_strand_damage_burden$Diagnosis <- factor(single_strand_damage_burden$Diagnosis,
                                                levels = c("Non-smoker_Control", "Smoker_Control", "Smoker_ASCVD"),
                                                labels = c("Non-smoker w/o ASCVD", "Smoker w/o ASCVD", "Smoker w/ ASCVD"))

# Remove cells with too small callable regions (bottom 5%)
single_strand_damage_burden <- single_strand_damage_burden[single_strand_damage_burden$callable_region_size_filter,]

# Mixed effect model
ssSNV_lm <- glmmTMB(
  snv ~ Age + Smoking*Depth + ASCVD +
    offset(log(ssSNV_effective_callable_region_size)) +
    (1 | Case_ID),
  data = single_strand_damage_burden,
  family = nbinom2
)
summary(ssSNV_lm) -> ssSNV_lm_summary
p_smoking_ssSNV <- ssSNV_lm_summary$coefficients$cond["SmokingSmoker", "Pr(>|z|)"]
p_ASCVD_ssSNV <- ssSNV_lm_summary$coefficients$cond["ASCVDASCVD", "Pr(>|z|)"]

pred <- ggpredict(ssSNV_lm, terms = "Age", 
                  condition = c(ssSNV_effective_callable_region_size = 1), 
                  type = "fixed")
rate_age_lm <- lm(predicted ~ x, data = pred)
ssSNV_age_slope <- coef(rate_age_lm)["x"]
ssSNV_age_intercept <- coef(rate_age_lm)["(Intercept)"]

# Plot
# SNV
ggplot(single_strand_damage_burden, 
       aes(x = Age, y = ssSNV_rate, color = Diagnosis)) +
  geom_point(size=0.5,shape=16) +
  geom_abline(intercept=ssSNV_age_intercept,slope=ssSNV_age_slope,
              size=0.5,color="lightskyblue") + 
  xlab("Age (yrs)") + ylab("ssDNA SNVs per bp") + 
  annotate(geom="text", x=0, y=2.85e-6*c(1, 0.925), 
           label = c(
             paste0("Smoking:~italic(P)==", scientific_10(p_smoking_ssSNV)),
             paste0("ASCVD:~italic(P)==", scientific_10(p_ASCVD_ssSNV))
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
    plot.title = element_text(color = "black", size = font_size, family = font),
    axis.line.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.title.y.left = element_text(margin = margin(r = -8))
  ) +
  scale_y_break(c(1.4e-6, 2.5e-6), scales = 0.25) +
  scale_y_continuous(labels = scientific_10, breaks = c(0, 5e-7, 1e-6, 2.5e-6, 3e-6), limits = c(0, 3e-6))
#ggsave("./figures/ssDNA_mut_manuscript_figures/figure4/ssSNV_burden.pdf", width = 5.4, height = 5, units = "cm")

# ------------------------------------------------------------------------------
# dsDNA burden
EC_raw=read.delim("./data/mutation_data/heart_EC.summary.mutation_rate.txt", header=T)
EC_filtered=EC_raw[!(grepl("down-sampling",EC_raw$Note)|grepl("redundant",EC_raw$Note)|grepl("excluded",EC_raw$Note)),]
EC_final_dsSNV=data.frame(EC_filtered$Cell_ID,EC_filtered$SCAN2_rate,as.factor(EC_filtered$Case_ID),EC_filtered$Age,EC_filtered$Gender,EC_filtered$Tissue_origin,EC_filtered$Clinical,EC_filtered$MAPD_HC,EC_filtered$CoV_HC)
colnames(EC_final_dsSNV)=c("Cell_ID","Somatic_rate","Case_ID","Age","Gender","Tissue_origin","Clinical","MAPD","CoV")
EC_final_dsSNV$sSNV_number=EC_final_dsSNV$Somatic_rate*genomic_factor
all_EC_cells=EC_final_dsSNV$Cell_ID[!is.na(EC_final_dsSNV$Somatic_rate)]

EC_final_dsSNV=EC_final_dsSNV[!is.na(EC_final_dsSNV$Somatic_rate),]
EC_final_dsSNV <- EC_final_dsSNV[EC_final_dsSNV$Clinical != "ASCVD",]

EC_final_dsSNV$Clinical=factor(EC_final_dsSNV$Clinical,levels=c("none","smoking","both"))
EC_final_dsSNV$Smoking=factor(EC_final_dsSNV$Clinical,
                              levels=c("none","smoking","both"),
                              labels=c("Non-smoker", "Smoker", "Smoker"))
EC_final_dsSNV$ASCVD=factor(EC_final_dsSNV$Clinical,
                            levels=c("none","smoking","both"),
                            labels=c("Control", "Control", "ASCVD"))
mutation_rates <- EC_final_dsSNV$Somatic_rate
EC_final_dsSNV$Diagnosis <- str_c(EC_final_dsSNV$Smoking, "_", EC_final_dsSNV$ASCVD)
EC_final_dsSNV$Diagnosis <- factor(EC_final_dsSNV$Diagnosis,
                                   levels = c("Non-smoker_Control", "Smoker_Control", "Smoker_ASCVD"),
                                   labels = c("Non-smoker w/o ASCVD", "Smoker w/o ASCVD", "Smoker w/ ASCVD"))


# mixed effect model
dsSNV_lm <- lmer(sSNV_number ~ Age + Smoking + ASCVD + (1|Case_ID),
                 data = EC_final_dsSNV,
                 REML=FALSE
)
summary(dsSNV_lm) -> dsSNV_lm_summary
p_smoking_dsSNV <- dsSNV_lm_summary$coefficients["SmokingSmoker", "Pr(>|t|)"]
p_ASCVD_dsSNV <- dsSNV_lm_summary$coefficients["ASCVDASCVD", "Pr(>|t|)"]
dsSNV_age_slope=dsSNV_lm_summary$coefficients["Age", "Estimate"]
dsSNV_age_intercept=dsSNV_lm_summary$coefficients["(Intercept)", "Estimate"]

# plot
ggplot(EC_final_dsSNV, 
       aes(x = Age, y = sSNV_number/genomic_factor/1e9, color = Diagnosis)) +
  geom_point(size=0.5,shape=16) +
  geom_abline(intercept=dsSNV_age_intercept/genomic_factor/1e9,
              slope=dsSNV_age_slope/genomic_factor/1e9,
              size=0.5,color="lightskyblue") + 
  xlab("Age (yrs)") + ylab("dsDNA SNVs per bp") + 
  annotate(geom="text", x=0, y=3.3e-6*c(1, 0.925), 
           label = c(
             paste0("Smoking:~italic(P)==", scientific_10(p_smoking_dsSNV)),
             paste0("ASCVD:~italic(P)==", scientific_10(p_ASCVD_dsSNV))
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
  scale_y_continuous(labels = scientific_10, n.breaks = 3, breaks = c(0, 1.5e-6, 3e-6))
#ggsave("./figures/ssDNA_mut_manuscript_figures/figure4/dsSNV_burden.pdf", width = 5, height = 4.5, units = "cm")

# ------------------------------------------------------------------------------
# Single-stranded damage spectrum (CVEC)
single_strand_damage_mat <- readRDS("./results/ssDNA_mut_spectrum/EC_single_strand_damage_mat_age_matched.rds")

# SNV
single_strand_SNV_mat_condition <- single_strand_damage_mat[["single_strand_SNV_mat_condition"]]
single_strand_SNV_mat_condition_diff <- single_strand_damage_mat[["single_strand_SNV_mat_condition_diff"]]
single_strand_SNV_mat_condition <- cbind(single_strand_SNV_mat_condition, single_strand_SNV_mat_condition_diff)
colnames(single_strand_SNV_mat_condition) <- 
  c("Non-smoker w/o\nASCVD (Age-matched)", "Smoker w/o ASCVD", "Smoker w/ ASCVD", 
    "Smoker - Non-smoker\n(w/o ASCVD)", "ASCVD - Non-ASCVD\n(Smoker)")
single_strand_SNV_mat_condition[single_strand_SNV_mat_condition < 0] <- 0

my.plot_96_profile_absolute(single_strand_SNV_mat_condition, 
                            condensed = T, facet_grid_scales = "free_x") + 
  theme_classic() +
  theme(
    strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    axis.text.x = element_text(color = "black", size = 2, family = font, angle = 90, vjust = 0.5),
    axis.ticks.x=element_blank(), 
    panel.spacing.x = unit(0.1, "lines"),
    strip.text.y = element_text(color = "black", size = 5, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
  ) +
  labs(y = "Absolute contribution (ssDNA SNVs per bp)", x = "Context") +
  scale_y_continuous(labels = scientific_10, n.breaks = 3)
ggsave("./figures/ssDNA_mut_manuscript_figures/figure4/ssSNV_spectrum.pdf", width = 9, height = 9, units = "cm")


# ------------------------------------------------------------------------------
# Single-stranded signature analysis
# COSMIC signature 
snv_signatures <- readRDS("../ref/COSMIC_signatures/COSMIC_v3.2_SBS_GRCh37.rds")
artefact_SBS_signatures <- c("SBS27", "SBS43", "SBS45", "SBS46", "SBS47", "SBS48", 
                             "SBS49", "SBS50", "SBS51", "SBS52", "SBS53", "SBS54", 
                             "SBS55", "SBS56", "SBS57", "SBS58", "SBS59", "SBS60", 
                             "SBS95")
snv_signatures <- snv_signatures[, !colnames(snv_signatures) %in% artefact_SBS_signatures]

# spectrum by cell
single_strand_damage_mat_cell <- readRDS("./results/ssDNA_mut_spectrum/EC_single_strand_damage_mat_cell.rds")
single_strand_SNV_mat_cell <- single_strand_damage_mat_cell$single_strand_SNV_mat_cell[, single_strand_damage_burden$Cell_ID]
single_strand_Indel_mat_cell <- single_strand_damage_mat_cell$single_strand_Indel_mat_cell[, single_strand_damage_burden$Cell_ID]

# SBS signature selection (use the smoking samples)
SNV_signature_fit <- fit_signatures_until_stable(single_strand_SNV_mat_condition, 
                                                 snv_signatures, max_delta = 0.01)
selected_signatures <- 
  c(SNV_signature_fit$per_sample$`Non-smoker w/o\nASCVD (Age-matched)`$selected_signatures,
    SNV_signature_fit$per_sample$`Smoker w/o ASCVD`$selected_signatures,
    SNV_signature_fit$per_sample$`Smoker w/ ASCVD`$selected_signatures)
selected_signatures <- unique(selected_signatures)

# signature fitting
SNV_signature_fit <- fit_to_signatures(single_strand_SNV_mat_cell, snv_signatures[,selected_signatures])
apply(SNV_signature_fit$contribution, 2, function(X){X/(sum(X)+1e-31)}) -> SNV_signature_fit
SNV_signature_fit <- data.frame(t(SNV_signature_fit))
SNV_signature_fit[,"Cell_ID"] <- row.names(SNV_signature_fit)
SNV_signature_fit <- right_join(single_strand_damage_burden, SNV_signature_fit)

# SBS4
SNV_signature_fit$SBS4_tot_count <- SNV_signature_fit$SBS4 * SNV_signature_fit$snv
SNV_signature_fit$SBS4_rate <- SNV_signature_fit$SBS4_tot_count / SNV_signature_fit$ssSNV_effective_callable_region_size

# mixed effect model
ssSBS4_lm <- glmmTMB(
  SBS4_tot_count ~ Age + Smoking*Depth + ASCVD +
    offset(log(ssSNV_effective_callable_region_size)) +
    (1 | Case_ID),
  data = SNV_signature_fit,
  family = nbinom2
)
summary(ssSBS4_lm) -> ssSBS4_lm_summary
p_smoking_ssSBS4 <- ssSBS4_lm_summary$coefficients$cond["SmokingSmoker", "Pr(>|z|)"]
p_ASCVD_ssSBS4 <- ssSBS4_lm_summary$coefficients$cond["ASCVDASCVD", "Pr(>|z|)"]

pred <- ggpredict(ssSBS4_lm, terms = "Age", 
                  condition = c(ssSNV_effective_callable_region_size = 1), 
                  type = "fixed")
rate_age_lm <- lm(predicted ~ x, data = pred)
ssSBS4_age_slope <- coef(rate_age_lm)["x"]
ssSBS4_age_intercept <- coef(rate_age_lm)["(Intercept)"]

# plot
ggplot(SNV_signature_fit,aes(x=Age,y=SBS4_rate,color=Diagnosis)) +
  geom_point(size=0.5,shape=16) +
  geom_abline(intercept=ssSBS4_age_intercept,slope=ssSBS4_age_slope,
              size=0.5,color="lightskyblue") + 
  xlab("Age (yrs)") + ylab("ssDNA SBS4 contribution per bp") + 
  annotate(geom="text", x=0, y=10e-7*c(1, 0.925), 
           label = c(
             paste0("Smoking:~italic(P)==", scientific_10(p_smoking_ssSBS4)),
             paste0("ASCVD:~italic(P)==", scientific_10(p_ASCVD_ssSBS4))
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
  scale_y_continuous(labels = scientific_10)
ggsave("./figures/ssDNA_mut_manuscript_figures/figure4/ssSBS4_burden.pdf", width = 5, height = 4.5, units = "cm")
