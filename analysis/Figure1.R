SSDH_REPO <- getwd()   # repository root; run these scripts from there
setwd("/n/data1/bch/genetics/lee/shulin/heart_endo")

.libPaths("/n/data1/bch/genetics/lee/shulin/heart_endo/Rpackages_O2_upgraded/")
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(scales)
library(reshape2)

.libPaths("/home/shm283/shulin/softwares/anaconda3/envs/mutationalpatterns/lib/R/library")
library(MutationalPatterns)

source(file.path(SSDH_REPO, "analysis/utils.R"))

# font
font = "Helvetica"
font_size = 7
axis_font_size = 6

# Custom function to format scientific notation with superscripts and handle zero
scientific_10 <- function(x) {
  if(x<=0.05){
    formatted_x <- sprintf("%.2e", x)  # Convert to scientific notation
    gsub("e\\+?", " %*% 10^", formatted_x)  # Replace 'e+' with ' * 10^'
  }else{
    sprintf("%.2f", x)
  }
}

# ------------------------------------------------------------------------------
# Single-stranded damage spectrum (control neurons + ECs)
neuron_single_strand_damage_mat <- readRDS("./results/ssDNA_mut_spectrum/neuron_single_strand_damage_mat.rds")
EC_single_strand_damage_mat <- readRDS("./results/ssDNA_mut_spectrum/EC_single_strand_damage_mat.rds")

# COSMIC signature 
snv_signatures <- readRDS("../ref/COSMIC_signatures/COSMIC_v3.2_SBS_GRCh37.rds")

# HiDEF-seq results
HiDEFseq_neuron_ssDNA_spectrum <- readRDS("./results/ssDNA_mut_spectrum/HiDEFseq_neuron_ssDNA_SNV_spectrum.rds")

# SNV
neuron_single_strand_SNV_mat <- neuron_single_strand_damage_mat[["single_strand_SNV_mat_condition"]]
EC_single_strand_SNV_mat <- EC_single_strand_damage_mat[["single_strand_SNV_mat_condition"]]
single_strand_SNV_mat_condition <- cbind(neuron_single_strand_SNV_mat[,"Control"], 
                                         EC_single_strand_SNV_mat[,"Control (Non-smoker)"], 
                                         snv_signatures[, "SBS30"])
single_strand_SNV_mat_condition <- apply(single_strand_SNV_mat_condition, 2, function(X){X/sum(X)})
colnames(single_strand_SNV_mat_condition) <- c("Neuron", "CVEC", "SBS30")

cos_sim_single_strand_SNV_sbs30 <- round(cos_sim_matrix(single_strand_SNV_mat_condition, snv_signatures[,"SBS30", drop=F])[1:2, 1], 2)
cos_sim_single_strand_SNV_HiDEFseq <- round(cos_sim_matrix(single_strand_SNV_mat_condition, HiDEFseq_neuron_ssDNA_spectrum)["Neuron", "HiDEF-seq"], 2)


my.plot_96_profile_absolute(single_strand_SNV_mat_condition, condensed = T) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = 4, family = font, angle = 90, vjust = 0.5),
    axis.ticks.x=element_blank(), 
    panel.spacing.x = unit(0.1, "lines"),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
  ) +
  labs(y = "Relative contribution of ssDNA SNV", x = "Context", 
       caption = paste0("Cosine similiaity to SBS30: ", cos_sim_single_strand_SNV_sbs30["Neuron"],
                        "\nto HiDEF-seq neuron ssDNA SNV: ", cos_sim_single_strand_SNV_HiDEFseq,
                        "\nCosine similiaity to SBS30: ", cos_sim_single_strand_SNV_sbs30["CVEC"]))
ggsave("./figures/ssDNA_mut_manuscript_figures/figure1/sSNV_spectrum.pdf", width = 12.5, height = 8, units = "cm")


# Indel
neuron_single_strand_Indel_mat <- neuron_single_strand_damage_mat[["single_strand_Indel_mat_condition"]]
EC_single_strand_Indel_mat <- EC_single_strand_damage_mat[["single_strand_Indel_mat_condition"]]
single_strand_Indel_mat_condition <- cbind(neuron_single_strand_Indel_mat[,"Control"], 
                                           EC_single_strand_Indel_mat[,"Control (Non-smoker)"])
single_strand_Indel_mat_condition <- apply(single_strand_Indel_mat_condition, 2, function(X){X/sum(X)})
colnames(single_strand_Indel_mat_condition) <- c("Neuron", "CVEC")

my.plot_indel_contexts(single_strand_Indel_mat_condition, condensed = T) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = 4, family = font, angle = 90, vjust = 0.5),
    axis.ticks.x = element_blank(), 
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    panel.spacing.x = unit(0.1, "lines"),
    legend.position = "none"
  ) +
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  labs(y = "Relative contribution of ssDNA Indel", x = "Context")
ggsave("./figures/ssDNA_mut_manuscript_figures/figure1/sIndel_spectrum.pdf", width = 12.5, height = 5.6, units = "cm")

# ------------------------------------------------------------------------------
# VAF of single-stranded damage
single_strand_mut_control_neuron <- readRDS("./results/ssDNA_mut/control_neuron_single_strand_damage.rds")
single_strand_mut_EC <- readRDS("./results/ssDNA_mut/EC_single_strand_damage.rds")

# read in metadata
single_strand_damage_burden <- readRDS("./results/ssDNA_mut_burden/EC_single_strand_damage_burden.rds")
control_EC_cell_id <- single_strand_damage_burden$Cell_ID[
  single_strand_damage_burden$Smoking == "Non-smoker" &
    single_strand_damage_burden$ASCVD == "Control"]
single_strand_mut_control_EC <- single_strand_mut_EC[single_strand_mut_EC$cell_id %in% control_EC_cell_id,]

single_strand_mut_control_neuron$cell_type <- "Neuron"
single_strand_mut_control_EC$cell_type <- "CVEC"
single_strand_mut <- rbind(single_strand_mut_control_neuron, single_strand_mut_control_EC)

# Calculate the corrected vaf
allele_1 <- abs(single_strand_mut$af - 0.5*single_strand_mut$ab)
allele_2 <- abs(single_strand_mut$af - 0.5*(1-single_strand_mut$ab))
ifelse(allele_1 < allele_2, 
       single_strand_mut$af / (single_strand_mut$ab) / 2,
       single_strand_mut$af / ((1-single_strand_mut$ab)) / 2) -> corrected_af
single_strand_mut$corrected_af <- corrected_af

# plot
single_strand_mut$muttype <- factor(single_strand_mut$muttype,
                                    levels = c("snv", "indel"), labels = c("SNV", "Indel"))
single_strand_mut$cell_type <- factor(single_strand_mut$cell_type,
                                    levels = c("Neuron", "CVEC"))

ggplot(single_strand_mut, aes(x = corrected_af, fill = muttype)) +
  geom_histogram(position = position_identity(), binwidth = 0.005) +
  scale_fill_manual(values = c("#fed339", "#7099e1")) +
  facet_wrap(.~cell_type, ncol=1, strip.position = "right") +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.key.size = unit(0.3, "line")
  )  +
  labs(x = "ssDNA damage VAF\n(corrected by allele balance)", y = "ssDNA damage count", fill = NULL) +
  scale_x_continuous(limits = c(0, 1))
ggsave("./figures/ssDNA_mut_manuscript_figures/figure1/corrected_VAF.pdf", width = 5, height = 6, units = "cm")


# ------------------------------------------------------------------------------
# Mis-classified rates of double-strand mutations  (control neurons)
neuron_misclassified_germline_mutation <- readRDS("./results/ssDNA_mut_error_estimation/Control_neuron_misclassified_rate.rds")
neuron_misclassified_germline_mutation$snv[,"muttype"] <- "snv"
colnames(neuron_misclassified_germline_mutation$snv)[27] <- "misclassified_dsDNA_calls_in_ss_calls"
neuron_misclassified_germline_mutation$indel[,"muttype"] <- "indel"
colnames(neuron_misclassified_germline_mutation$indel)[27] <- "misclassified_dsDNA_calls_in_ss_calls"
neuron_misclassified_germline_mutation <- rbind(neuron_misclassified_germline_mutation$snv[, c("muttype", "misclassified_dsDNA_calls_in_ss_calls")],
                                                neuron_misclassified_germline_mutation$indel[, c("muttype", "misclassified_dsDNA_calls_in_ss_calls")])
neuron_misclassified_germline_mutation$cell_type <- "Neuron"

EC_misclassified_germline_mutation <- readRDS("./results/ssDNA_mut_error_estimation/EC_misclassified_rate.rds")
EC_misclassified_germline_mutation$snv[,"muttype"] <- "snv"
colnames(EC_misclassified_germline_mutation$snv)[31] <- "misclassified_dsDNA_calls_in_ss_calls"
EC_misclassified_germline_mutation$indel[,"muttype"] <- "indel"
colnames(EC_misclassified_germline_mutation$indel)[31] <- "misclassified_dsDNA_calls_in_ss_calls"
EC_misclassified_germline_mutation <- rbind(EC_misclassified_germline_mutation$snv[, c("muttype", "misclassified_dsDNA_calls_in_ss_calls")],
                                            EC_misclassified_germline_mutation$indel[, c("muttype", "misclassified_dsDNA_calls_in_ss_calls")])
EC_misclassified_germline_mutation$cell_type <- "CVEC"

misclassified_germline_mutation <- rbind(neuron_misclassified_germline_mutation, 
                                         EC_misclassified_germline_mutation)
misclassified_germline_mutation$muttype <- factor(misclassified_germline_mutation$muttype, 
                                                         levels = c("snv", "indel"), labels = c("SNV", "Indel"))
misclassified_germline_mutation$cell_type <- factor(misclassified_germline_mutation$cell_type,
                                      levels = c("Neuron", "CVEC"))


ggplot(misclassified_germline_mutation, aes(x = muttype, y = misclassified_dsDNA_calls_in_ss_calls, fill = muttype)) +
  geom_violin(linewidth = 0.3) +
  geom_jitter(size = 0.05) + 
  facet_wrap(.~cell_type, ncol=1, strip.position = "right") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#fed339", "#7099e1")) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )  +
  labs(x = NULL, y = "Estimated rate of dsDNA mutation\nmisclassified as ssDNA damage in cells")
ggsave("./figures/ssDNA_mut_manuscript_figures/figure1/misclassification_rate.pdf", width = 5, height = 5, units = "cm")

