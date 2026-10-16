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
library(rstatix)
library(ggh4x)
library(ggbreak)
library(lme4)
library(lmerTest)
library(ggeffects)

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

indel_signatures <- readRDS("../ref/COSMIC_signatures/COSMIC_v3.2_ID_GRCh37.rds")
indel_signature_name <- readRDS("../ref/COSMIC_signatures/COSMIC_ID_type_name.rds")

genomic_factor <- 5.845


table(single_strand_damage_burden$Diagnosis)
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
ssIndel_lm <- glmmTMB(
  indel ~ Age + Smoking*Depth + ASCVD +
    offset(log(ssIndel_effective_callable_region_size)) +
    (1 | Case_ID),
  data = single_strand_damage_burden,
  family = nbinom2
)
summary(ssIndel_lm) -> ssIndel_lm_summary
p_smoking_ssIndel <- ssIndel_lm_summary$coefficients$cond["SmokingSmoker", "Pr(>|z|)"]
p_ASCVD_ssIndel <- ssIndel_lm_summary$coefficients$cond["ASCVDASCVD", "Pr(>|z|)"]

pred <- ggpredict(ssIndel_lm, terms = "Age", 
                  condition = c(ssIndel_effective_callable_region_size = 1), 
                  type = "fixed")
rate_age_lm <- lm(predicted ~ x, data = pred)
ssIndel_age_slope <- coef(rate_age_lm)["x"]
ssIndel_age_intercept <- coef(rate_age_lm)["(Intercept)"]

# Plot
# Indel
ggplot(single_strand_damage_burden, 
       aes(x = Age, y = ssIndel_rate, color = Diagnosis)) +
  geom_point(size=0.5,shape=16) +
  geom_abline(intercept=ssIndel_age_intercept,slope=ssIndel_age_slope,
              size=0.5,color="lightskyblue") + 
  xlab("Age (yrs)") + ylab("ssDNA Indels per bp") + 
  annotate(geom="text", x=0, y=1.1e-6*c(1, 0.925), 
           label = c(
             paste0("Smoking:~italic(P)==", scientific_10(p_smoking_ssIndel)),
             paste0("ASCVD:~italic(P)==", scientific_10(p_ASCVD_ssIndel))
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
#ggsave("./figures/ssDNA_mut_manuscript_figures/figure3/ssIndel_burden.pdf", width = 5, height = 4.5, units = "cm")

# ------------------------------------------------------------------------------
# Single-stranded damage spectrum (CVEC)
single_strand_damage_mat <- readRDS("./results/ssDNA_mut_spectrum/EC_single_strand_damage_mat_age_matched_ASCVD_Nonsmoker.rds")

# Indel
single_strand_Indel_mat_condition <- single_strand_damage_mat[["single_strand_Indel_mat_condition"]]
single_strand_Indel_mat_condition_diff <- single_strand_damage_mat[["single_strand_Indel_mat_condition_diff"]]
single_strand_Indel_mat_condition <- cbind(single_strand_Indel_mat_condition, single_strand_Indel_mat_condition_diff)
colnames(single_strand_Indel_mat_condition) <- 
  c("Non-smoker w/o ASCVD (age-matched)", "Non-smoker w/ ASCVD", "Smoker w/o ASCVD", "Smoker w/ ASCVD",
    "ASCVD - Non-ASCVD (Non-smoker)", "Smoker - Non-smoker (w/o ASCVD)", "ASCVD - Non-ASCVD (Smoker)")
single_strand_Indel_mat_condition <- single_strand_Indel_mat_condition[,
                                                                       c("Non-smoker w/o ASCVD (age-matched)", "Non-smoker w/ ASCVD", "Smoker w/o ASCVD", "Smoker w/ ASCVD",
                                                                         "Smoker - Non-smoker (w/o ASCVD)", "ASCVD - Non-ASCVD (Non-smoker)", "ASCVD - Non-ASCVD (Smoker)")
]

# Signature analysis
# COSMIC signature
single_strand_Indel_mat_condition_signature <- fit_to_signatures(single_strand_Indel_mat_condition, indel_signatures)
single_strand_Indel_mat_condition_signature_contribution <- melt(single_strand_Indel_mat_condition_signature$contribution)
colnames(single_strand_Indel_mat_condition_signature_contribution) <- c("signature", "diagnosis", "contribution")
single_strand_Indel_mat_condition_signature_contribution$diagnosis <- factor(
  single_strand_Indel_mat_condition_signature_contribution$diagnosis,
  levels = rev(c("Non-smoker w/o ASCVD (age-matched)", "Non-smoker w/ ASCVD", "Smoker w/o ASCVD", "Smoker w/ ASCVD",
                 "Smoker - Non-smoker (w/o ASCVD)", "ASCVD - Non-ASCVD (Non-smoker)", "ASCVD - Non-ASCVD (Smoker)"))
)

ggplot(single_strand_Indel_mat_condition_signature_contribution, aes(x = diagnosis, y = signature, fill = contribution)) + 
  geom_tile(color = "black", linewidth = 0.2) + # Add white borders for cells
  scale_fill_gradient2(low = "white", high = "red", # Custom diverging color scale
                       name = "Signature\ncontribution\nper bp",
                       labels = scientific_10) +
  labs(title = "ssDNA Indel", x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 90, vjust = 0.5, hjust = 1),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.title = element_text(color = "black", size = axis_font_size, family = font),
    legend.key.size = unit(0.5, "line"),
    plot.title = element_text(color = "black", size = font_size, family = font),
    axis.title = element_blank(), # Remove axis titles as they are redundant with labels
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    legend.position = "right"
  ) +
  coord_flip() +
  scale_x_discrete(position = "top")
ggsave("./figures/ssDNA_mut_manuscript_figures/figure3/ssIndel_COSMIC_signature_decomposition.pdf", width = 10, height = 4.3, units = "cm")

# Spectrum
single_strand_Indel_mat_condition[single_strand_Indel_mat_condition < 0] <- 0
my.plot_indel_contexts(single_strand_Indel_mat_condition, 
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
  scale_y_continuous(labels = scientific_10)
ggsave("./figures/ssDNA_mut_manuscript_figures/figure3/ssIndel_spectrum.pdf", width = 9, height = 12, units = "cm")

my.plot_indel_contexts(single_strand_Indel_mat_condition, 
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
  labs(y = "Absolute contribution (ssDNA Indels per bp)", x = "Context") +
  scale_y_continuous(labels = scientific_10)
ggsave("./figures/ssDNA_mut_manuscript_figures/figure3/ssIndel_spectrum_diff_nonsmoker_ASCVD.pdf", width = 9, height = 12, units = "cm")

# ------------------------------------------------------------------------------
# Single-stranded signature analysis
# spectrum by cell
single_strand_damage_mat_cell <- readRDS("./results/ssDNA_mut_spectrum/EC_single_strand_damage_mat_cell.rds")
single_strand_SNV_mat_cell <- single_strand_damage_mat_cell$single_strand_SNV_mat_cell[, single_strand_damage_burden$Cell_ID]
single_strand_Indel_mat_cell <- single_strand_damage_mat_cell$single_strand_Indel_mat_cell[, single_strand_damage_burden$Cell_ID]

# Signature fitting
Indel_signature_fit <- fit_to_signatures(single_strand_Indel_mat_cell, indel_signatures)
apply(Indel_signature_fit$contribution, 2, function(X){X/(sum(X)+1e-31)}) -> Indel_signature_fit
Indel_signature_fit <- data.frame(t(Indel_signature_fit))
Indel_signature_fit[,"Cell_ID"] <- row.names(Indel_signature_fit)
Indel_signature_fit <- left_join(single_strand_damage_burden, Indel_signature_fit)

# ID4
Indel_signature_fit$ID4_tot_count <- Indel_signature_fit$ID4 * Indel_signature_fit$indel
Indel_signature_fit$ID4_rate <- Indel_signature_fit$ID4_tot_count / Indel_signature_fit$ssIndel_effective_callable_region_size

# mixed effect model
ssID4_lm <- glmmTMB(
  ID4_tot_count ~ Age + Smoking*Depth + ASCVD +
    offset(log(ssIndel_effective_callable_region_size)) +
    (1 | Case_ID),
  data = Indel_signature_fit,
  family = nbinom2
)
summary(ssID4_lm) -> ssID4_lm_summary
p_smoking_ssID4 <- ssID4_lm_summary$coefficients$cond["SmokingSmoker", "Pr(>|z|)"]
p_ASCVD_ssID4 <- ssID4_lm_summary$coefficients$cond["ASCVDASCVD", "Pr(>|z|)"]

pred <- ggpredict(ssID4_lm, terms = "Age", 
                  condition = c(ssIndel_effective_callable_region_size = 1), 
                  type = "fixed")
rate_age_lm <- lm(predicted ~ x, data = pred)
ssID4_age_slope <- coef(rate_age_lm)["x"]
ssID4_age_intercept <- coef(rate_age_lm)["(Intercept)"]

# plot
ggplot(Indel_signature_fit,aes(x=Age,y=ID4_rate,color=Diagnosis)) +
  geom_point(size=0.5,shape=16) +
  geom_abline(intercept=ssID4_age_slope,slope=ssID4_age_slope,
              size=0.5,color="lightskyblue") + 
  xlab("Age (yrs)") + ylab("ssDNA ID4 contribution per bp") + 
  annotate(geom="text", x=0, y=4.95e-7*c(1, 0.965), 
           label = c(
             paste0("Smoking:~italic(P)==", scientific_10(p_smoking_ssID4)),
             paste0("ASCVD:~italic(P)==", scientific_10(p_ASCVD_ssID4))
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
  scale_y_break(c(2e-7, 4.5e-7), scales = 0.25) +
  scale_y_continuous(labels = scientific_10, breaks = c(0, 1e-7, 2e-7, 4.5e-7, 5e-7), limits = c(0, 5.05e-7))
ggsave("./figures/ssDNA_mut_manuscript_figures/figure3/ssID4_burden.pdf", width = 5.4, height = 5, units = "cm")

# ------------------------------------------------------------------------------
# Double-stranded damage burden (CVEC)
# Indel
EC_raw=read.delim("./data/mutation_data/heart_EC.summary.mutation_rate.txt", header=T)
EC_filtered=EC_raw[!(grepl("down-sampling",EC_raw$Note)|grepl("redundant",EC_raw$Note)|grepl("excluded",EC_raw$Note)),]
EC_final_dsIndel=data.frame(EC_filtered$Cell_ID,EC_filtered$SCAN2_indel_rate,as.factor(EC_filtered$Case_ID),EC_filtered$Age,EC_filtered$Gender,EC_filtered$Tissue_origin,EC_filtered$Clinical,EC_filtered$MAPD_HC,EC_filtered$CoV_HC)
colnames(EC_final_dsIndel)=c("Cell_ID","Somatic_rate","Case_ID","Age","Gender","Tissue_origin","Clinical","MAPD","CoV")
EC_final_dsIndel$sSNV_number=EC_final_dsIndel$Somatic_rate*genomic_factor
EC_final_dsIndel$sIndel_number=EC_final_dsIndel$Somatic_rate*genomic_factor
EC_cell_list_dsIndel=EC_final_dsIndel$Cell_ID[!is.na(EC_final_dsIndel$Somatic_rate)]

EC_final_dsIndel=EC_final_dsIndel[!is.na(EC_final_dsIndel$Somatic_rate),]

EC_final_dsIndel <- EC_final_dsIndel[EC_final_dsIndel$Clinical != "ASCVD",]
EC_final_dsIndel$Clinical=factor(EC_final_dsIndel$Clinical,levels=c("none","smoking","both"))
EC_final_dsIndel$Smoking=factor(EC_final_dsIndel$Clinical,
                                levels=c("none","smoking","both"),
                                labels=c("Non-smoker", "Smoker", "Smoker"))
EC_final_dsIndel$ASCVD=factor(EC_final_dsIndel$Clinical,
                              levels=c("none","smoking","both"),
                              labels=c("Control", "Control", "ASCVD"))
EC_final_dsIndel$Diagnosis <- str_c(EC_final_dsIndel$Smoking, "_", EC_final_dsIndel$ASCVD)
EC_final_dsIndel$Diagnosis <- factor(EC_final_dsIndel$Diagnosis,
                                     levels = c("Non-smoker_Control", "Smoker_Control", "Smoker_ASCVD"),
                                     labels = c("Non-smoker w/o ASCVD", "Smoker w/o ASCVD", "Smoker w/ ASCVD"))
indel_mutation_rates <- EC_final_dsIndel$Somatic_rate

# mixed effect model
dsIndel_lm <- lmer(sIndel_number ~ Age + Smoking + ASCVD + (1|Case_ID),
                   data = EC_final_dsIndel,
                   REML=FALSE
)
summary(dsIndel_lm) -> dsIndel_lm_summary
p_smoking_dsIndel <- dsIndel_lm_summary$coefficients["SmokingSmoker", "Pr(>|t|)"]
p_ASCVD_dsIndel <- dsIndel_lm_summary$coefficients["ASCVDASCVD", "Pr(>|t|)"]
dsIndel_age_slope=dsIndel_lm_summary$coefficients["Age", "Estimate"]
dsIndel_age_intercept=dsIndel_lm_summary$coefficients["(Intercept)", "Estimate"]

# plot
ggplot(EC_final_dsIndel, 
       aes(x = Age, y = sIndel_number/genomic_factor/1e9, color = Diagnosis)) +
  geom_point(size=0.5,shape=16) +
  geom_abline(intercept=dsIndel_age_intercept/genomic_factor/1e9,
              slope=dsIndel_age_slope/genomic_factor/1e9,
              size=0.5,color="lightskyblue") + 
  xlab("Age (yrs)") + ylab("dsDNA Indels per bp") + 
  annotate(geom="text", x=0, y=2.2e-7*c(1, 0.925), 
           label = c(
             paste0("Smoking:~italic(P)==", scientific_10(p_smoking_dsIndel)),
             paste0("ASCVD:~italic(P)==", scientific_10(p_ASCVD_dsIndel))
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
#ggsave("./figures/ssDNA_mut_manuscript_figures/figure3/dsIndel_burden.pdf", width = 5, height = 4.5, units = "cm")

# ------------------------------------------------------------------------------
# Double-stranded damage spectrum (CVEC)
# Indel
# read in vcf files of the spectrum
vcf_files <- sprintf("./data/mutation_data/vcf_spectrum/%s.vcf", EC_cell_list_dsIndel)
vcfs_indel <- read_vcfs_as_granges(vcf_files, EC_cell_list_dsIndel, ref_genome, type="indel")
indel <- get_indel_context(vcfs_indel, ref_genome)
mut_mat_indel <- count_indel_contexts(indel)

# some cell with 0 sindel called
matrix(0, nrow = dim(mut_mat_indel)[1], ncol = sum(!EC_final_dsIndel$Cell_ID %in% colnames(mut_mat_indel))) -> indel_counts_0_cell
colnames(indel_counts_0_cell) <- EC_final_dsIndel$Cell_ID[!EC_final_dsIndel$Cell_ID %in% colnames(mut_mat_indel)]
rownames(indel_counts_0_cell) <- row.names(mut_mat_indel)
cbind(mut_mat_indel, indel_counts_0_cell) -> mut_mat_indel
mut_mat_indel <- mut_mat_indel[, EC_final_dsIndel$Cell_ID]

absolute_mut_mat_indel=t(t(mut_mat_indel)/(colSums(mut_mat_indel)+1e-31)*indel_mutation_rates*genomic_factor)
absolute_mut_mat_indel_smoking <- cbind(
  rowMeans(absolute_mut_mat_indel[,EC_final_dsIndel$Cell_ID[EC_final_dsIndel$Diagnosis == "Non-smoker w/o ASCVD"]]),
  rowMeans(absolute_mut_mat_indel[,EC_final_dsIndel$Cell_ID[EC_final_dsIndel$Diagnosis == "Non-smoker w/ ASCVD"]]),
  rowMeans(absolute_mut_mat_indel[,EC_final_dsIndel$Cell_ID[EC_final_dsIndel$Diagnosis == "Smoker w/o ASCVD"]]),
  rowMeans(absolute_mut_mat_indel[,EC_final_dsIndel$Cell_ID[EC_final_dsIndel$Diagnosis == "Smoker w/ ASCVD"]]))
colnames(absolute_mut_mat_indel_smoking) <- c("Non-smoker w/o ASCVD", "Non-smoker w/ ASCVD", "Smoker w/o ASCVD", "Smoker w/ ASCVD")

difference <- cbind(
  absolute_mut_mat_indel_smoking[,colnames(absolute_mut_mat_indel_smoking)=="Smoker w/o ASCVD"]-
    absolute_mut_mat_indel_smoking[,colnames(absolute_mut_mat_indel_smoking)=="Non-smoker w/o ASCVD"],
  absolute_mut_mat_indel_smoking[,colnames(absolute_mut_mat_indel_smoking)=="Non-smoker w/ ASCVD"]-
    absolute_mut_mat_indel_smoking[,colnames(absolute_mut_mat_indel_smoking)=="Non-smoker w/o ASCVD"],
  absolute_mut_mat_indel_smoking[,colnames(absolute_mut_mat_indel_smoking)=="Smoker w/ ASCVD"]-
    absolute_mut_mat_indel_smoking[,colnames(absolute_mut_mat_indel_smoking)=="Smoker w/o ASCVD"]
)
absolute_mut_mat_indel_smoking <- cbind(absolute_mut_mat_indel_smoking, difference)
colnames(absolute_mut_mat_indel_smoking) <-
  c("Non-smoker w/o\nASCVD (Age-matched)", "Non-smoker w/ ASCVD", "Smoker w/o ASCVD", "Smoker w/ ASCVD", 
    "Smoker - Non-smoker\n(w/o ASCVD)", "ASCVD - Non-ASCVD\n(Non-smoker)", "ASCVD - Non-ASCVD\n(Smoker)")

absolute_mut_mat_indel_smoking[absolute_mut_mat_indel_smoking < 0] <- 0
absolute_mut_mat_indel_smoking <- absolute_mut_mat_indel_smoking/genomic_factor/1e9

# plot
my.plot_indel_contexts(absolute_mut_mat_indel_smoking[, c("ASCVD - Non-ASCVD\n(Non-smoker)", 
                                                          "ASCVD - Non-ASCVD\n(Smoker)"), drop=F], 
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
  labs(y = "Absolute contribution (dsDNA Indels per bp)", x = "Context") +
  scale_y_continuous(labels = scientific_10)
ggsave("./figures/ssDNA_mut_manuscript_figures/figure3/dsIndel_spectrum.pdf", width = 9.2, height = 4.5, units = "cm")

# ID4
fit_Indel_signature <- fit_to_signatures(absolute_mut_mat_indel,indel_signatures)
fit_Indel_signature <- data.frame(t(fit_Indel_signature$contribution))
fit_Indel_signature[,"Cell_ID"] <- row.names(fit_Indel_signature)
fit_Indel_signature <- dplyr::right_join(EC_final_dsIndel, fit_Indel_signature)

# mixed effect model
dsID4_lm <- lmer(ID4 ~ Age + Smoking + ASCVD + (1|Case_ID),
                 data = fit_Indel_signature,
                 REML = FALSE
)
summary(dsID4_lm) -> dsID4_lm_summary
p_smoking_dsID4 <- dsID4_lm_summary$coefficients["SmokingSmoker", "Pr(>|t|)"]
p_ASCVD_dsID4 <- dsID4_lm_summary$coefficients["ASCVDASCVD", "Pr(>|t|)"]
dsID4_age_slope=dsID4_lm_summary$coefficients["Age", "Estimate"]
dsID4_age_intercept=dsID4_lm_summary$coefficients["(Intercept)", "Estimate"]

ggplot(fit_Indel_signature, 
       aes(x = Age, y = ID4/genomic_factor/1e9, color = Diagnosis)) +
  geom_point(size=0.5,shape=16) +
  geom_abline(intercept=dsID4_age_intercept/genomic_factor/1e9,
              slope=dsID4_age_slope/genomic_factor/1e9,
              size=0.5,color="lightskyblue") + 
  xlab("Age (yrs)") + ylab("dsDNA ID4 contribution per bp") + 
  annotate(geom="text", x=0, y=7.5e-8*c(1, 0.925), 
           label = c(
             paste0("Smoking:~italic(P)==", scientific_10(p_smoking_dsID4)),
             paste0("ASCVD:~italic(P)==", scientific_10(p_ASCVD_dsID4))
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
ggsave("./figures/ssDNA_mut_manuscript_figures/figure3/dsID4_burden.pdf", width = 4.8, height = 4.5, units = "cm")

# ------------------------------------------------------------------------------
# ID3 and ID4 spectrum

scientific_10_2 <- function(x) {
  sapply(x, function(v) {
    if (is.na(v)) {
      NA_character_
    } else if (v == 0) {
      "0"
    } else {
      parse(text = gsub("e", " %*% 10^", sprintf("%.2e", v)))
    }
  })
}


ID3_4 <- indel_signatures[, c("ID3", "ID4")]
indel_signature_name$scan2_name == rownames(ID3_4)
rownames(ID3_4) <- indel_signature_name$mutationalpattern_name
my.plot_indel_contexts(ID3_4, 
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
  labs(y = "Relative\ncontribution", x = "Context") +
  scale_y_continuous(labels = scientific_10_2)
ggsave("./figures/ssDNA_mut_manuscript_figures/figure3/ID3_4_spectrum.pdf", width = 9.2, height = 4.5, units = "cm")
