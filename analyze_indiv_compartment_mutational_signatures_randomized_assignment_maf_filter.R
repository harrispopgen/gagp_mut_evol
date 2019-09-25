# Analyze mutational signatures per individual (lots of PCAs)

setwd("~/Documents/Harris_project_directory")
# setwd("/net/harris/vol1/project/primate_ervs")

options(stringsAsFactors = F)
library(ggplot2)
library(ggridges)
library(plyr)
library(RSvgDevice)
library(ggfortify)
library(scales)
library(grid)
library(gridExtra)
library(cowplot)


source("./common.R")

# old_control_trinuc_composition <- read.table("./hg18_nonrepetitive_neutral_trinucs.txt", header = T)
control_trinuc_composition <- read.table("./hg18_control/hg18_control_3mer_content.txt", header = T)
early_rep_trinuc_composition <- read.table("./hg18_replication_timing_q3_liftovers/Homo_hg18_replication_timing_q3_3mer_content.txt", header = T)
late_rep_trinuc_composition <- read.table("./hg18_replication_timing_q0_liftovers/Homo_hg18_replication_timing_q0_3mer_content.txt", header = T)
erv_trinuc_composition <- read.table("./hg18_all_erv_trinucs.txt", header = T)
heterochromatin_trinuc_composition <- read.table("./chromHMM_tracks/hg18_nonrepetitive_chromHMM_heterochromatin_trinuc_composition.txt", header = T)
line_trinuc_composition <- read.table("./hg18_all_line_liftovers/Homo_hg18_all_line_3mer_content.txt", header = T)
early_rep_nonrepetitive_trinuc_composition <- read.table("./hg18_replication_timing_q3_nonrepetitive_liftovers/Homo_hg18_replication_timing_q3_nonrepetitive_3mer_content.txt", header = T)
early_rep_repetitive_trinuc_composition <- read.table("./hg18_replication_timing_q3_repetitive_liftovers/Homo_hg18_replication_timing_q3_repetitive_3mer_content.txt", header = T)
late_rep_nonrepetitive_trinuc_composition <- read.table("./hg18_replication_timing_q0_nonrepetitive_liftovers/Homo_hg18_replication_timing_q0_nonrepetitive_3mer_content.txt", header = T)
late_rep_repetitive_trinuc_composition <- read.table("./hg18_replication_timing_q0_repetitive_liftovers/Homo_hg18_replication_timing_q0_repetitive_3mer_content.txt", header = T)
low_recom_trinuc_composition <- read.table("./decode_recombination_hotspots_compartments/recombination_compartment_1_of_10_3mer_content.txt", header = T)
high_recom_trinuc_composition <- read.table("./decode_recombination_hotspots_compartments/recombination_compartment_10_of_10_3mer_content.txt", header = T)
chimp_low_recom_trinuc_composition <- read.table("./panMap_recom_hotspots/chimp_recom_coldspots_hg18_3mer_content.txt", header = T)
chimp_high_recom_trinuc_composition <- read.table("./panMap_recom_hotspots/chimp_recom_hotspots_hg18_3mer_content.txt", header = T)
# pericentromere_trinuc_composition <- read.table("./pericentromeric_regions/hg18_pericentromeric_regions_3mer_content.txt", header = T)
# subtelomere_trinuc_composition <- read.table("./subtelomeric_regions/hg18_subtelomeric_regions_3mer_content.txt", header = T)
erv_hmc_high_trinuc_composition <- read.table("./hydroxymethylation_compartment/hg18_all_erv_hmc_high_3mer_content.txt", header = T)
erv_hmc_low_trinuc_composition <- read.table("./hydroxymethylation_compartment/hg18_all_erv_hmc_low_3mer_content.txt", header = T)
maternal_hotspots_trinuc_composition <- read.table("./maternal_hotspots/hg18_maternal_hotspots_zero_indexed_mb_3mer_content.txt", header = T)
# maternal_hotspots_chr16_trinuc_composition <- read.table("./maternal_hotspots/hg18_maternal_hotspots_chr16_3mer_content.txt", header = T)
cpg_islands_trinuc_composition <- read.table("./cpg_islands/hg18_cpgIsland_3mer_content.txt", header = T)

control_spectra <- load_indiv_3mer_spectrum("hg18_control", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
erv_spectra <- load_indiv_3mer_spectrum("hg18_all_erv", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
early_rep_spectra <- load_indiv_3mer_spectrum("hg18_replication_timing_q3", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
late_rep_spectra <- load_indiv_3mer_spectrum("hg18_replication_timing_q0", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
heterochromatin_spectra <- load_indiv_3mer_spectrum("hg18_nonrepetitive_chromHMM_heterochromatin", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
line_spectra <- load_indiv_3mer_spectrum("hg18_all_line", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
early_rep_repetitive_spectra <- load_indiv_3mer_spectrum("hg18_replication_timing_q3_repetitive", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
early_rep_nonrepetitive_spectra <- load_indiv_3mer_spectrum("hg18_replication_timing_q3_nonrepetitive", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
late_rep_repetitive_spectra <- load_indiv_3mer_spectrum("hg18_replication_timing_q0_repetitive", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
late_rep_nonrepetitive_spectra <- load_indiv_3mer_spectrum("hg18_replication_timing_q0_nonrepetitive", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
# low_recom_spectra <- load_indiv_3mer_spectrum("recombination_compartment_1_of_10", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
# high_recom_spectra <- load_indiv_3mer_spectrum("recombination_compartment_10_of_10", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
# chimp_low_recom_spectra <- load_indiv_3mer_spectrum("chimp_recom_coldspots_hg18", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
# chimp_high_recom_spectra <- load_indiv_3mer_spectrum("chimp_recom_hotspots_hg18", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
# pericentromere_spectra <- load_indiv_3mer_spectrum("hg18_pericentromeric_regions")
# subtelomere_spectra <- load_indiv_3mer_spectrum("hg18_subtelomeric_regions")
maternal_hotspots_spectra <- load_indiv_3mer_spectrum("hg18_maternal_hotspots_zero_indexed_mb", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
# maternal_hotspots_chr16_spectra <- load_indiv_3mer_spectrum("hg18_maternal_hotspots_chr16", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
cpg_islands_spectra <- load_indiv_3mer_spectrum("hg18_cpgIsland", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")

control_spectra_reweight <- reweight_indiv_spectra(control_spectra, control_trinuc_composition, control_trinuc_composition)
erv_spectra_reweight <- reweight_indiv_spectra(erv_spectra, erv_trinuc_composition, control_trinuc_composition)
early_rep_spectra_reweight <- reweight_indiv_spectra(early_rep_spectra, early_rep_trinuc_composition, control_trinuc_composition)
late_rep_spectra_reweight <- reweight_indiv_spectra(late_rep_spectra, late_rep_trinuc_composition, control_trinuc_composition)
heterochromatin_spectra_reweight <- reweight_indiv_spectra(heterochromatin_spectra, heterochromatin_trinuc_composition, control_trinuc_composition)
line_spectra_reweight <- reweight_indiv_spectra(line_spectra, line_trinuc_composition, control_trinuc_composition)
early_rep_repetitive_spectra_reweight <- reweight_indiv_spectra(early_rep_repetitive_spectra, early_rep_repetitive_trinuc_composition, control_trinuc_composition)
early_rep_nonrepetitive_spectra_reweight <- reweight_indiv_spectra(early_rep_nonrepetitive_spectra, early_rep_nonrepetitive_trinuc_composition, control_trinuc_composition)
late_rep_repetitive_spectra_reweight <- reweight_indiv_spectra(late_rep_repetitive_spectra, late_rep_repetitive_trinuc_composition, control_trinuc_composition)
late_rep_nonrepetitive_spectra_reweight <- reweight_indiv_spectra(late_rep_nonrepetitive_spectra, late_rep_nonrepetitive_trinuc_composition, control_trinuc_composition)
# low_recom_spectra_reweight <- reweight_indiv_spectra(low_recom_spectra, low_recom_trinuc_composition, control_trinuc_composition)
# high_recom_spectra_reweight <- reweight_indiv_spectra(high_recom_spectra, high_recom_trinuc_composition, control_trinuc_composition)
# chimp_low_recom_spectra_reweight <- reweight_indiv_spectra(chimp_low_recom_spectra, chimp_low_recom_trinuc_composition, control_trinuc_composition)
# chimp_high_recom_spectra_reweight <- reweight_indiv_spectra(chimp_high_recom_spectra, chimp_high_recom_trinuc_composition, control_trinuc_composition)
# pericentromere_spectra_reweight <- reweight_indiv_spectra(pericentromere_spectra, pericentromere_trinuc_composition, control_trinuc_composition)
# subtelomere_spectra_reweight <- reweight_indiv_spectra(subtelomere_spectra, subtelomere_trinuc_composition, control_trinuc_composition)
maternal_hotspots_spectra_reweight <- reweight_indiv_spectra(maternal_hotspots_spectra, maternal_hotspots_trinuc_composition, control_trinuc_composition)
# maternal_hotspots_chr16_spectra_reweight <- reweight_indiv_spectra(maternal_hotspots_chr16_spectra, maternal_hotspots_chr16_trinuc_composition, control_trinuc_composition)
cpg_islands_spectra_reweight <- reweight_indiv_spectra(cpg_islands_spectra, cpg_islands_trinuc_composition, control_trinuc_composition)

spectra_by_species_repetitive_split_reweight <-
  lapply(
    species,
    function(s)
      rbind(
        subset(control_spectra_reweight, startsWith(rownames(control_spectra_reweight), s)),
        subset(erv_spectra_reweight, startsWith(rownames(erv_spectra_reweight), s)),
        subset(early_rep_repetitive_spectra_reweight, startsWith(rownames(early_rep_repetitive_spectra_reweight), s)),
        subset(early_rep_nonrepetitive_spectra_reweight, startsWith(rownames(early_rep_nonrepetitive_spectra_reweight), s)),
        subset(late_rep_repetitive_spectra_reweight, startsWith(rownames(late_rep_repetitive_spectra_reweight), s)),
        subset(late_rep_nonrepetitive_spectra_reweight, startsWith(rownames(late_rep_nonrepetitive_spectra_reweight), s)),
        subset(heterochromatin_spectra_reweight, startsWith(rownames(heterochromatin_spectra_reweight), s)),
        subset(line_spectra_reweight, startsWith(rownames(line_spectra_reweight), s))
        )
  )
names(spectra_by_species_repetitive_split_reweight) <- species

indiv_df_by_species_repetitive_split <-
  lapply(
    species,
    function(s){
      n_indiv = nrow(subset(control_spectra, startsWith(rownames(control_spectra), s)))
      return(
        data.frame(
          "indiv" = rownames(spectra_by_species_repetitive_split_reweight[[s]]),
          "subspecies" = sub("-.*", "", rownames(spectra_by_species_repetitive_split_reweight[[s]])),
          "spectrum" = c(rep("control", n_indiv),
                         rep("erv", n_indiv),
                         rep("early_rep_repetitive", n_indiv),
                         rep("early_rep_nonrepetitive", n_indiv),
                         rep("late_rep_repetitive", n_indiv),
                         rep("late_rep_nonrepetitive", n_indiv),
                         rep("heterochromatin", n_indiv),
                         rep("line", n_indiv)
                         ))
      )
    }
  )
names(indiv_df_by_species_repetitive_split) <- species

# eigen_windows(
#   t(
#     as.matrix(
#       spectra_by_species_repetitive_split_reweight$Homo
#       )
#     ), k=3, win = 96
#   )

pca_by_species_repetitive_split_reweight <-
  lapply(
    spectra_by_species_repetitive_split_reweight,
    function(spectra)
      prcomp(
        spectra,
        scale. = T
      )
  )
names(pca_by_species_repetitive_split_reweight) <- species

ggsave(
  "./pca_by_species/all_spectrum_repetitive_split_pca_reweight_maf_filter_Homo_pc1_pc2_20190827.pdf",
  ggplot(
    cbind(pca_by_species_repetitive_split_reweight$Homo$x, indiv_df_by_species_repetitive_split$Homo),
    aes(x=PC1, y=PC2, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[4]) +
    scale_fill_manual(values=subspecies_colors[4]) +
    xlab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Homo, 1)) +
    ylab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Homo, 2)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/all_spectrum_repetitive_split_pca_reweight_maf_filter_Homo_pc1_pc3_20190827.pdf",
  ggplot(
    cbind(pca_by_species_repetitive_split_reweight$Homo$x, indiv_df_by_species_repetitive_split$Homo),
    aes(x=PC1, y=PC3, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[4]) +
    scale_fill_manual(values=subspecies_colors[4]) +
    xlab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Homo, 1)) +
    ylab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Homo, 3)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/all_spectrum_repetitive_split_pca_reweight_maf_filter_negative_PC1_Homo_pc1_pc2_20190831.pdf",
  ggplot(
    cbind(pca_by_species_repetitive_split_reweight$Homo$x, indiv_df_by_species_repetitive_split$Homo),
    aes(x=(-1*PC1), y=PC2, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[4]) +
    scale_fill_manual(values=subspecies_colors[4]) +
    xlab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Homo, 1)) +
    ylab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Homo, 2)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/all_spectrum_repetitive_split_pca_reweight_maf_filter_negative_PC1_Homo_pc1_pc3_20190831.pdf",
  ggplot(
    cbind(pca_by_species_repetitive_split_reweight$Homo$x, indiv_df_by_species_repetitive_split$Homo),
    aes(x=(-1*PC1), y=PC3, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[4]) +
    scale_fill_manual(values=subspecies_colors[4]) +
    xlab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Homo, 1)) +
    ylab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Homo, 3)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/all_spectrum_repetitive_split_pca_reweight_maf_filter_Pan_troglodytes_pc1_pc2_20190827.pdf",
  ggplot(
    cbind(pca_by_species_repetitive_split_reweight$Pan_troglodytes$x, indiv_df_by_species_repetitive_split$Pan_troglodytes),
    aes(x=PC1, y=PC2, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[c(6:9)]) +
    scale_fill_manual(values=subspecies_colors[c(6:9)]) +
    xlab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pan_troglodytes, 1)) +
    ylab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pan_troglodytes, 2)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/all_spectrum_repetitive_split_pca_reweight_maf_filter_Pan_troglodytes_pc1_pc3_20190827.pdf",
  ggplot(
    cbind(pca_by_species_repetitive_split_reweight$Pan_troglodytes$x, indiv_df_by_species_repetitive_split$Pan_troglodytes),
    aes(x=PC1, y=PC3, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[c(6:9)]) +
    scale_fill_manual(values=subspecies_colors[c(6:9)]) +
    xlab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pan_troglodytes, 1)) +
    ylab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pan_troglodytes, 3)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/all_spectrum_repetitive_split_pca_reweight_maf_filter_Gorilla_pc1_pc2_20190827.pdf",
  ggplot(
    cbind(pca_by_species_repetitive_split_reweight$Gorilla$x, indiv_df_by_species_repetitive_split$Gorilla),
    aes(x=PC1, y=PC2, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[c(1:3)]) +
    scale_fill_manual(values=subspecies_colors[c(1:3)]) +
    xlab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Gorilla, 1)) +
    ylab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Gorilla, 2)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/all_spectrum_repetitive_split_pca_reweight_maf_filter_Gorilla_pc1_pc3_20190827.pdf",
  ggplot(
    cbind(pca_by_species_repetitive_split_reweight$Gorilla$x, indiv_df_by_species_repetitive_split$Gorilla),
    aes(x=PC1, y=PC3, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[c(1:3)]) +
    scale_fill_manual(values=subspecies_colors[c(1:3)]) +
    xlab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Gorilla, 1)) +
    ylab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Gorilla, 3)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/all_spectrum_repetitive_split_pca_reweight_maf_filter_Pan_paniscus_pc1_pc2_20190827.pdf",
  ggplot(
    cbind(pca_by_species_repetitive_split_reweight$Pan_paniscus$x, indiv_df_by_species_repetitive_split$Pan_paniscus),
    aes(x=PC1, y=PC2, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[c(5)]) +
    scale_fill_manual(values=subspecies_colors[c(5)]) +
    xlab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pan_paniscus, 1)) +
    ylab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pan_paniscus, 2)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/all_spectrum_repetitive_split_pca_reweight_maf_filter_Pan_paniscus_pc1_pc3_20190827.pdf",
  ggplot(
    cbind(pca_by_species_repetitive_split_reweight$Pan_paniscus$x, indiv_df_by_species_repetitive_split$Pan_paniscus),
    aes(x=PC1, y=PC3, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[c(5)]) +
    scale_fill_manual(values=subspecies_colors[c(5)]) +
    xlab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pan_paniscus, 1)) +
    ylab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pan_paniscus, 3)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/all_spectrum_repetitive_split_pca_reweight_maf_filter_Pongo_abelii_pc1_pc2_20190827.pdf",
  ggplot(
    cbind(pca_by_species_repetitive_split_reweight$Pongo_abelii$x, indiv_df_by_species_repetitive_split$Pongo_abelii),
    aes(x=PC1, y=PC2, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[c(10)]) +
    scale_fill_manual(values=subspecies_colors[c(10)]) +
    xlab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pongo_abelii, 1)) +
    ylab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pongo_abelii, 2)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/all_spectrum_repetitive_split_pca_reweight_maf_filter_negative_PC2_Pongo_abelii_pc1_pc2_20190831.pdf",
  ggplot(
    cbind(pca_by_species_repetitive_split_reweight$Pongo_abelii$x, indiv_df_by_species_repetitive_split$Pongo_abelii),
    aes(x=PC1, y=(-1*PC2), group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[c(10)]) +
    scale_fill_manual(values=subspecies_colors[c(10)]) +
    xlab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pongo_abelii, 1)) +
    ylab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pongo_abelii, 2)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/all_spectrum_repetitive_split_pca_reweight_maf_filter_Pongo_abelii_pc1_pc3_20190827.pdf",
  ggplot(
    cbind(pca_by_species_repetitive_split_reweight$Pongo_abelii$x, indiv_df_by_species_repetitive_split$Pongo_abelii),
    aes(x=PC1, y=PC3, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[c(10)]) +
    scale_fill_manual(values=subspecies_colors[c(10)]) +
    xlab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pongo_abelii, 1)) +
    ylab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pongo_abelii, 3)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/all_spectrum_repetitive_split_pca_reweight_maf_filter_negative_PC3_Pongo_abelii_pc1_pc3_20190831.pdf",
  ggplot(
    cbind(pca_by_species_repetitive_split_reweight$Pongo_abelii$x, indiv_df_by_species_repetitive_split$Pongo_abelii),
    aes(x=PC1, y=(-1*PC3), group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[c(10)]) +
    scale_fill_manual(values=subspecies_colors[c(10)]) +
    xlab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pongo_abelii, 1)) +
    ylab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pongo_abelii, 3)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/all_spectrum_repetitive_split_pca_reweight_maf_filter_Pongo_pygmaeus_pc1_pc2_20190827.pdf",
  ggplot(
    cbind(pca_by_species_repetitive_split_reweight$Pongo_pygmaeus$x, indiv_df_by_species_repetitive_split$Pongo_pygmaeus),
    aes(x=PC1, y=PC2, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[c(11)]) +
    scale_fill_manual(values=subspecies_colors[c(11)]) +
    xlab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pongo_pygmaeus, 1)) +
    ylab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pongo_pygmaeus, 2)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/all_spectrum_repetitive_split_pca_reweight_maf_filter_Pongo_pygmaeus_pc1_pc3_20190827.pdf",
  ggplot(
    cbind(pca_by_species_repetitive_split_reweight$Pongo_pygmaeus$x, indiv_df_by_species_repetitive_split$Pongo_pygmaeus),
    aes(x=PC1, y=PC3, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[c(11)]) +
    scale_fill_manual(values=subspecies_colors[c(11)]) +
    xlab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pongo_pygmaeus, 1)) +
    ylab(generate_pc_axis_label(pca_by_species_repetitive_split_reweight$Pongo_pygmaeus, 3)),
  width = 7, height = 4)

ggsave("./ridgeplot_all_control_maf_filter_20190827.pdf",
       plot_indiv_compartment_chisq(control_spectra,
                                    control_trinuc_composition,
                                    control_spectra,
                                    control_trinuc_composition,
                                    "NRNC Chi2 Distribution", species)
)
ggsave("./ridgeplot_subset_control_maf_filter_20190827.pdf",
       plot_indiv_compartment_chisq(control_spectra,
                                    control_trinuc_composition,
                                    control_spectra,
                                    control_trinuc_composition,
                                    "NRNC Chi2 Distribution",
                                    c("Homo", "Pan_troglodytes", "Gorilla"))
)

indiv_df <-
  data.frame(
    "indiv" = rownames(control_spectra_reweight),
    "subspecies" = sub("-.*", "", rownames(control_spectra_reweight))
  )
indiv_df$species <-
  sapply(indiv_df$subspecies,
         function(x)
           if(startsWith(x, "Gorilla")) "Gorilla" else 
             if(startsWith(x, "Pan_troglodytes")) "Pan_troglodytes" else x)
control_pca <-
  prcomp(
    control_spectra_reweight,
    scale. = T
  )

indiv_df_control_rep_time <- 
  data.frame(
    "indiv" = 
      c(rownames(control_spectra_reweight),
        rownames(early_rep_spectra_reweight),
        rownames(late_rep_spectra_reweight)),
    "subspecies"=
      c(
        sub("-.*", "", rownames(control_spectra_reweight)),
        sub("-.*", "", rownames(early_rep_spectra_reweight)),
        sub("-.*", "", rownames(late_rep_spectra_reweight))
      ),
    "spectrum" = c(rep("control", nrow(control_spectra_reweight)),
                   rep("early_rep", nrow(early_rep_spectra_reweight)),
                   rep("late_rep", nrow(late_rep_spectra_reweight)))
  )
indiv_df_control_rep_time$species <-
  sapply(indiv_df_control_rep_time$subspecies,
         function(x)
           if(startsWith(x, "Gorilla")) "Gorilla" else 
             if(startsWith(x, "Pan_troglodytes")) "Pan_troglodytes" else x)
control_early_late_rep_pca <-
  prcomp(
    rbind(
      control_spectra_reweight,
      early_rep_spectra_reweight,
      late_rep_spectra_reweight
    ),
    scale. = T
  )

spectra_by_species_replication_timing <-
  lapply(
    species,
    function(s)
      rbind(
        subset(control_spectra_reweight, startsWith(rownames(control_spectra_reweight), s)),
        subset(early_rep_spectra_reweight, startsWith(rownames(control_spectra_reweight), s)),
        subset(late_rep_spectra_reweight, startsWith(rownames(control_spectra_reweight), s))
      )
  )
names(spectra_by_species_replication_timing) <- species

control_early_late_rep_pca_by_species <-
  lapply(
    species,
    function(s)
      prcomp(
        spectra_by_species_replication_timing[[s]],
        scale. = T
      )
  )
names(control_early_late_rep_pca_by_species) <- species

indiv_df_by_species_control_early_late_rep <-
  lapply(
    species,
    function(s){
      n_indiv = nrow(subset(control_spectra, startsWith(rownames(control_spectra), s)))
      return(
        data.frame(
          "indiv" = rownames(spectra_by_species_replication_timing[[s]]),
          "subspecies" = sub("-.*", "", rownames(spectra_by_species_replication_timing[[s]])),
          "spectrum" = c(rep("NCNR", n_indiv),
                         rep("early_rep", n_indiv),
                         rep("late_rep", n_indiv)
          ))
      )
    }
  )
names(indiv_df_by_species_control_early_late_rep) <- species

ggsave(
  "./pca_by_species/rep_timing_pca_reweight_maf_filter_Gorilla_pc1_pc2_20190827.pdf",
  ggplot(
    cbind(control_early_late_rep_pca_by_species$Gorilla$x, indiv_df_by_species_control_early_late_rep$Gorilla),
    aes(x=PC1, y=PC2, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(17, 15, 16))+
    scale_color_manual(values=subspecies_colors[c(1:3)]) +
    scale_fill_manual(values=subspecies_colors[c(1:3)]) +
    xlab(generate_pc_axis_label(control_early_late_rep_pca_by_species$Gorilla, 1)) +
    ylab(generate_pc_axis_label(control_early_late_rep_pca_by_species$Gorilla, 2)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/rep_timing_pca_reweight_maf_filter_Homo_pc1_pc2_20190827.pdf",
  ggplot(
    cbind(control_early_late_rep_pca_by_species$Homo$x, indiv_df_by_species_control_early_late_rep$Homo),
    aes(x=PC1, y=PC2, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(17, 15, 16))+
    scale_color_manual(values=subspecies_colors[4]) +
    scale_fill_manual(values=subspecies_colors[4]) +
    xlab(generate_pc_axis_label(control_early_late_rep_pca_by_species$Homo, 1)) +
    ylab(generate_pc_axis_label(control_early_late_rep_pca_by_species$Homo, 2)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/rep_timing_pca_reweight_maf_filter_Pan_troglodytes_pc1_pc2_20190827.pdf",
  ggplot(
    cbind(control_early_late_rep_pca_by_species$Pan_troglodytes$x, indiv_df_by_species_control_early_late_rep$Pan_troglodytes),
    aes(x=PC1, y=PC2, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(17, 15, 16))+
    scale_color_manual(values=subspecies_colors[c(6:9)]) +
    scale_fill_manual(values=subspecies_colors[c(6:9)]) +
    xlab(generate_pc_axis_label(control_early_late_rep_pca_by_species$Pan_troglodytes, 1)) +
    ylab(generate_pc_axis_label(control_early_late_rep_pca_by_species$Pan_troglodytes, 2)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/rep_timing_pca_reweight_maf_filter_Pan_paniscus_pc1_pc2_20190913.pdf",
  ggplot(
    cbind(control_early_late_rep_pca_by_species$Pan_paniscus$x, indiv_df_by_species_control_early_late_rep$Pan_paniscus),
    aes(x=PC1, y=PC2, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(17, 15, 16))+
    scale_color_manual(values=subspecies_colors[5]) +
    scale_fill_manual(values=subspecies_colors[5]) +
    xlab(generate_pc_axis_label(control_early_late_rep_pca_by_species$Pan_paniscus, 1)) +
    ylab(generate_pc_axis_label(control_early_late_rep_pca_by_species$Pan_paniscus, 2)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/rep_timing_pca_reweight_maf_filter_Pongo_abelii_pc1_pc2_20190913.pdf",
  ggplot(
    cbind(control_early_late_rep_pca_by_species$Pongo_abelii$x, indiv_df_by_species_control_early_late_rep$Pongo_abelii),
    aes(x=PC1, y=PC2, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(17, 15, 16))+
    scale_color_manual(values=subspecies_colors[10]) +
    scale_fill_manual(values=subspecies_colors[10]) +
    xlab(generate_pc_axis_label(control_early_late_rep_pca_by_species$Pongo_abelii, 1)) +
    ylab(generate_pc_axis_label(control_early_late_rep_pca_by_species$Pongo_abelii, 2)),
  width = 7, height = 4)

ggsave(
  "./pca_by_species/rep_timing_pca_reweight_maf_filter_Pongo_pygmaeus_pc1_pc2_20190913.pdf",
  ggplot(
    cbind(control_early_late_rep_pca_by_species$Pongo_pygmaeus$x, indiv_df_by_species_control_early_late_rep$Pongo_pygmaeus),
    aes(x=PC1, y=PC2, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(17, 15, 16))+
    scale_color_manual(values=subspecies_colors[11]) +
    scale_fill_manual(values=subspecies_colors[11]) +
    xlab(generate_pc_axis_label(control_early_late_rep_pca_by_species$Pongo_pygmaeus, 1)) +
    ylab(generate_pc_axis_label(control_early_late_rep_pca_by_species$Pongo_pygmaeus, 2)),
  width = 7, height = 4)


# indiv_df_control_recom <- 
#   data.frame(
#     "indiv" = 
#       c(rownames(control_spectra_reweight),
#         rownames(low_recom_spectra_reweight),
#         rownames(high_recom_spectra_reweight)),
#     "subspecies"=
#       c(
#         sub("-.*", "", rownames(control_spectra_reweight)),
#         sub("-.*", "", rownames(low_recom_spectra_reweight)),
#         sub("-.*", "", rownames(high_recom_spectra_reweight))
#       ),
#     "spectrum" = c(rep("control", nrow(control_spectra_reweight)),
#                    rep("low_recom", nrow(low_recom_spectra_reweight)),
#                    rep("high_recom", nrow(high_recom_spectra_reweight)))
#   )

# indiv_df_control_recom$species <-
#   sapply(indiv_df_control_recom$subspecies,
#          function(x)
#            if(startsWith(x, "Gorilla")) "Gorilla" else 
#              if(startsWith(x, "Pan_troglodytes")) "Pan_troglodytes" else x)
# 
# control_recom_pca <-
#   prcomp(
#     rbind(
#       control_spectra_reweight,
#       low_recom_spectra_reweight,
#       high_recom_spectra_reweight
#     ),
#     scale. = T
#   )

# indiv_df_control_chimp_recom <-
#   data.frame(
#     "indiv" =
#       c(rownames(control_spectra_reweight),
#         rownames(chimp_low_recom_spectra_reweight),
#         rownames(chimp_high_recom_spectra_reweight)),
#     "subspecies"=
#       c(
#         sub("-.*", "", rownames(control_spectra_reweight)),
#         sub("-.*", "", rownames(chimp_low_recom_spectra_reweight)),
#         sub("-.*", "", rownames(chimp_high_recom_spectra_reweight))
#       ),
#     "spectrum" = c(rep("control", nrow(control_spectra_reweight)),
#                    rep("low_recom", nrow(chimp_low_recom_spectra_reweight)),
#                    rep("high_recom", nrow(chimp_high_recom_spectra_reweight)))
#   )
# 
# indiv_df_control_chimp_recom$species <-
#   sapply(indiv_df_control_chimp_recom$subspecies,
#          function(x)
#            if(startsWith(x, "Gorilla")) "Gorilla" else
#              if(startsWith(x, "Pan_troglodytes")) "Pan_troglodytes" else x)

# control_recom_chimp_pca <-
#   prcomp(
#     rbind(
#       control_spectra_reweight,
#       chimp_low_recom_spectra_reweight,
#       chimp_high_recom_spectra_reweight
#     ),
#     scale. = T
#   )

control_erv_pca <-
  prcomp(
    rbind(
      control_spectra_reweight,
      erv_spectra_reweight
    ),
    scale. = T
  )

indiv_df_control_erv <-
  data.frame(
    "indiv" =
      c(rownames(control_spectra_reweight),
        rownames(erv_spectra_reweight)),
    "subspecies"=
      c(
        sub("-.*", "", rownames(control_spectra_reweight)),
        sub("-.*", "", rownames(erv_spectra_reweight))
      ),
    "spectrum" = c(rep("control", nrow(control_spectra_reweight)),
                   rep("erv", nrow(erv_spectra_reweight)))
  )

indiv_df_control_erv$species <-
  sapply(indiv_df_control_erv$subspecies,
         function(x)
           if(startsWith(x, "Gorilla")) "Gorilla" else
             if(startsWith(x, "Pan_troglodytes")) "Pan_troglodytes" else x)


heterochromatin_erv_pca <-
  prcomp(
    rbind(
      heterochromatin_spectra_reweight,
      erv_spectra_reweight
    ),
    scale. = T
  )

indiv_df_heterochromatin_erv <-
  data.frame(
    "indiv" =
      c(rownames(heterochromatin_spectra_reweight),
        rownames(erv_spectra_reweight)),
    "subspecies"=
      c(
        sub("-.*", "", rownames(heterochromatin_spectra_reweight)),
        sub("-.*", "", rownames(erv_spectra_reweight))
      ),
    "spectrum" = c(rep("heterochromatin", nrow(heterochromatin_spectra_reweight)),
                   rep("erv", nrow(erv_spectra_reweight)))
  )

indiv_df_heterochromatin_erv$species <-
  sapply(indiv_df_heterochromatin_erv$subspecies,
         function(x)
           if(startsWith(x, "Gorilla")) "Gorilla" else
             if(startsWith(x, "Pan_troglodytes")) "Pan_troglodytes" else x)

maternal_hotspots_pca <-
  prcomp(
    maternal_hotspots_spectra_reweight,
    scale. = T
  )

# maternal_hotspots_chr16_pca <-
#   prcomp(
#     maternal_hotspots_chr16_spectra_reweight,
#     scale. = T
#   )

# control_maternal_hotspots_chr16_pca <-
#   prcomp(
#     rbind(
#       control_spectra_reweight,
#       maternal_hotspots_chr16_spectra_reweight
#     ),
#     scale. = T
#   )

control_maternal_hotspots_pca <-
  prcomp(
    rbind(
      control_spectra_reweight,
      maternal_hotspots_spectra_reweight
    ),
    scale. = T
  )

indiv_df_control_maternal_hotspots <-
  data.frame(
    "indiv" =
      c(rownames(control_spectra_reweight),
        rownames(maternal_hotspots_spectra_reweight)),
    "subspecies"=
      c(
        sub("-.*", "", rownames(control_spectra_reweight)),
        sub("-.*", "", rownames(maternal_hotspots_spectra_reweight))
      ),
    "spectrum" = c(rep("control", nrow(control_spectra_reweight)),
                   rep("maternal_hotspots", nrow(maternal_hotspots_spectra_reweight)))
  )

indiv_df_control_maternal_hotspots$species <-
  sapply(indiv_df_control_maternal_hotspots$subspecies,
         function(x)
           if(startsWith(x, "Gorilla")) "Gorilla" else
             if(startsWith(x, "Pan_troglodytes")) "Pan_troglodytes" else x)

# set.seed(1)
downsampled_indivs <-
  unlist(
    sapply(
      unique(indiv_df$species),
      function(x){
        sp_subset = subset(indiv_df$indiv, indiv_df$species == x)
        return(
          sample(
            subset(indiv_df$indiv, indiv_df$species == x),
            size = min(length(sp_subset), 9)
          )
        )
      }
    )
  )
# downsampled_control_pca <-
#   prcomp(
#     subset(
#       control_spectra_reweight,
#       rownames(control_spectra_reweight) %in% downsampled_indivs
#     ),
#     scale. = T
#   )
# indiv_df_downsampled <-
#   data.frame(
#     "indiv" = rownames(subset(control_spectra_reweight, rownames(control_spectra_reweight) %in% downsampled_indivs)),
#     "subspecies" = sub("-.*", "", rownames(subset(control_spectra_reweight, rownames(control_spectra_reweight) %in% downsampled_indivs)))
#   )
# indiv_df_downsampled$species <-
#   sapply(indiv_df_downsampled$subspecies,
#          function(x)
#            if(startsWith(x, "Gorilla")) "Gorilla" else 
#              if(startsWith(x, "Pan_troglodytes")) "Pan_troglodytes" else x)


control_cpg_islands_pca <-
  prcomp(
    rbind(
      control_spectra_reweight,
      cpg_islands_spectra_reweight
    ),
    scale. = T
  )

indiv_df_control_cpg_islands <-
  data.frame(
    "indiv" =
      c(rownames(control_spectra_reweight),
        rownames(cpg_islands_spectra_reweight)),
    "subspecies"=
      c(
        sub("-.*", "", rownames(control_spectra_reweight)),
        sub("-.*", "", rownames(cpg_islands_spectra_reweight))
      ),
    "spectrum" = c(rep("control", nrow(control_spectra_reweight)),
                   rep("cpg_islands", nrow(cpg_islands_spectra_reweight)))
  )

indiv_df_control_cpg_islands$species <-
  sapply(indiv_df_control_cpg_islands$subspecies,
         function(x)
           if(startsWith(x, "Gorilla")) "Gorilla" else
             if(startsWith(x, "Pan_troglodytes")) "Pan_troglodytes" else x)

# ggsave(
#   "./pca_all_species/pca_control_downsampled_pc1_pc2_20190827.pdf", 
#   ggplot(
#     cbind(downsampled_control_pca$x, indiv_df_downsampled), 
#     aes(x=PC1, y=PC2, group=species)) +
#     geom_point(aes(fill=species), color = "black", pch=21, size=5) +
#     scale_fill_manual(values=species_colors) +
#     xlab(generate_pc_axis_label(downsampled_control_pca, 1)) +
#     ylab(generate_pc_axis_label(downsampled_control_pca, 2)),
#   width=7, height=4
#   )
# 
# ggsave(
#   "./pca_all_species/pca_control_downsampled_pc1_pc3_20190827.pdf", 
#   ggplot(
#     cbind(downsampled_control_pca$x, indiv_df_downsampled), 
#     aes(x=PC1, y=PC3, group=species)) +
#     geom_point(aes(fill=species), color = "black", pch=21, size=5) +
#     scale_fill_manual(values=species_colors) +
#     xlab(generate_pc_axis_label(downsampled_control_pca, 1)) +
#     ylab(generate_pc_axis_label(downsampled_control_pca, 3)),
#   width=7, height=4
# )

ggsave("./pca_all_species/pca_control_subspecies_pc1_pc2_maf_filter_20190827.pdf",
       ggplot(
         cbind(control_pca$x, indiv_df),
         aes(x=PC1, y=PC2, group=subspecies)) +
         geom_point(aes(fill=subspecies), colour="black", pch=21, size=5) +
         scale_fill_manual(values=subspecies_colors) +
         xlab(generate_pc_axis_label(control_pca, 1)) +
         ylab(generate_pc_axis_label(control_pca, 2)),
       width=7, height=4)

# Donald is a hybrid w chimp with c chimp admixture! he's closer to c chimps in the PC space
# mean(
#   subset(
#     control_pca$x[,c("PC2")], 
#     startsWith(rownames(control_pca$x), "Pan_troglodytes_troglodytes")
#   )
# )
# 
# mean(
#   subset(
#     control_pca$x[,c("PC2")], 
#     startsWith(rownames(control_pca$x), "Pan_troglodytes_verus") &
#       !endsWith(rownames(control_pca$x), "Donald")
#   )
# )

ggsave(
  "./pca_all_species/pca_control_subspecies_pc1_pc3_maf_filter_20190827.pdf",
  ggplot(
    cbind(control_pca$x, indiv_df),
    aes(x=PC1, y=PC3, group=subspecies)) +
    geom_point(aes(fill=subspecies), colour="black", pch=21, size=5) +
    scale_fill_manual(values=subspecies_colors) +
    xlab(generate_pc_axis_label(control_pca, 1)) +
    ylab(generate_pc_axis_label(control_pca, 3)),
  width = 7, height = 4
)

# ggsave(
#   "./pca_all_species/pca_control_recom_subspecies_pc1_pc2_20190603.pdf", 
#   ggplot(
#     cbind(control_recom_pca$x, indiv_df_control_recom), 
#     aes(x=PC1, y=PC2, group=subspecies)) +
#     geom_point(aes(color=subspecies, shape=spectrum), size=5) +
#     scale_color_manual(values=subspecies_colors) +
#     xlab(generate_pc_axis_label(control_recom_pca, 1)) +
#     ylab(generate_pc_axis_label(control_recom_pca, 2)),
#   width = 7, height = 4)

# ggsave(
#   "./pca_all_species/pca_control_recom_subspecies_pc1_pc3_20190603.pdf", 
#   ggplot(
#     cbind(control_recom_pca$x, indiv_df_control_recom), 
#     aes(x=PC1, y=PC3, group=subspecies)) +
#     geom_point(aes(color=subspecies, shape=spectrum), size=5) +
#     scale_color_manual(values=subspecies_colors) +
#     xlab(generate_pc_axis_label(control_recom_pca, 1)) +
#     ylab(generate_pc_axis_label(control_recom_pca, 3)),
#   width = 7, height = 4)

# ggsave(
#   "./pca_all_species/pca_control_chimp_recom_subspecies_pc1_pc2_20190604.pdf",
#   ggplot(
#     cbind(control_recom_chimp_pca$x, indiv_df_control_chimp_recom),
#     aes(x=PC1, y=PC2, group=subspecies)) +
#     geom_point(aes(color=subspecies, shape=spectrum), size=5) +
#     scale_color_manual(values=subspecies_colors) +
#     xlab(generate_pc_axis_label(control_recom_chimp_pca, 1)) +
#     ylab(generate_pc_axis_label(control_recom_chimp_pca, 2)),
#   width = 7, height = 4)

# ggsave(
#   "./pca_all_species/pca_control_chimp_recom_subspecies_pc1_pc3_20190604.pdf",
#   ggplot(
#     cbind(control_recom_chimp_pca$x, indiv_df_control_chimp_recom),
#     aes(x=PC1, y=PC3, group=subspecies)) +
#     geom_point(aes(color=subspecies, shape=spectrum), size=5) +
#     scale_color_manual(values=subspecies_colors) +
#     xlab(generate_pc_axis_label(control_recom_chimp_pca, 1)) +
#     ylab(generate_pc_axis_label(control_recom_chimp_pca, 3)),
#   width = 7, height = 4)

# ggsave(
#   "./pca_all_species/pca_control_repressed_regions_pc1_pc2_20190509.pdf", 
#   ggplot(
#     cbind(control_subtelomere_pericentromere_pca$x, indiv_df_control_subtelomere_pericentromere),
#     aes(x=PC1, y=PC2, group=subspecies)) +
#     geom_point(aes(color=subspecies, shape=spectrum), size=5) +
#     scale_color_manual(values=subspecies_colors) +
#     xlab(generate_pc_axis_label(control_subtelomere_pericentromere_pca, 1)) +
#     ylab(generate_pc_axis_label(control_subtelomere_pericentromere_pca, 2)),
#   width = 7, height = 4)

# ggsave(
#   "./pca_all_species/pca_control_repressed_regions_pc1_pc3_20190509.pdf",
#   ggplot(
#     cbind(control_subtelomere_pericentromere_pca$x, indiv_df_control_subtelomere_pericentromere), 
#     aes(x=PC1, y=PC3, group=subspecies)) +
#     geom_point(aes(color=subspecies, shape=spectrum), size=5) +
#     scale_color_manual(values=subspecies_colors) +
#     xlab(generate_pc_axis_label(control_subtelomere_pericentromere_pca, 1)) +
#     ylab(generate_pc_axis_label(control_subtelomere_pericentromere_pca, 3)),
#   width = 7, height = 4)

ggsave(
  "./pca_all_species/pca_control_rep_timing_pc1_pc2_maf_filter_20190827.pdf",
  ggplot(
    cbind(control_early_late_rep_pca$x, indiv_df_control_rep_time), 
    aes(x=PC1, y=PC2, group=subspecies)) +
    geom_point(aes(color=subspecies, shape=spectrum), size=5) +
    scale_color_manual(values=subspecies_colors) +
    xlab(generate_pc_axis_label(control_early_late_rep_pca, 1)) +
    ylab(generate_pc_axis_label(control_early_late_rep_pca, 2)),
  width = 7, height = 4)

ggsave(
  "./pca_all_species/pca_control_rep_timing_pc1_pc3_maf_filter_20190827.pdf",
  ggplot(
    cbind(control_early_late_rep_pca$x, indiv_df_control_rep_time), 
    aes(x=PC1, y=PC3, group=subspecies)) +
    geom_point(aes(color=subspecies, shape=spectrum), size=5) +
    scale_color_manual(values=subspecies_colors) +
    xlab(generate_pc_axis_label(control_early_late_rep_pca, 1)) +
    ylab(generate_pc_axis_label(control_early_late_rep_pca, 3)),
  width = 7, height = 4)

# ggsave(
#   "./pca_all_species/pca_control_erv_subspecies_pc1_pc2_maf_filter_20190702.pdf", 
#   ggplot(
#     cbind(control_erv_pca$x, indiv_df_control_erv), 
#     aes(x=PC1, y=PC2, group=subspecies)) +
#     geom_point(aes(color=subspecies, shape=spectrum), size=5) +
#     scale_color_manual(values=subspecies_colors) +
#     xlab(generate_pc_axis_label(control_erv_pca, 1)) +
#     ylab(generate_pc_axis_label(control_erv_pca, 2)),
#   width = 7, height = 4)

# ggsave(
#   "./pca_all_species/pca_control_erv_subspecies_pc1_pc3_maf_filter_20190702.pdf", 
#   ggplot(
#     cbind(control_erv_pca$x, indiv_df_control_erv), 
#     aes(x=PC1, y=PC3, group=subspecies)) +
#     geom_point(aes(color=subspecies, shape=spectrum), size=5) +
#     scale_color_manual(values=subspecies_colors) +
#     xlab(generate_pc_axis_label(control_erv_pca, 1)) +
#     ylab(generate_pc_axis_label(control_erv_pca, 3)),
#   width = 7, height = 4)


ggsave(
  "./pca_all_species/pca_control_cpg_islands_subspecies_pc1_pc2_maf_filter_20190827.pdf", 
  ggplot(
    cbind(control_cpg_islands_pca$x, indiv_df_control_cpg_islands), 
    aes(x=PC1, y=PC2, group=subspecies)) +
    geom_point(aes(color=subspecies, shape=spectrum), size=5) +
    scale_color_manual(values=subspecies_colors) +
    xlab(generate_pc_axis_label(control_cpg_islands_pca, 1)) +
    ylab(generate_pc_axis_label(control_cpg_islands_pca, 2)),
  width = 7, height = 4)

# ggsave(
#   "./pca_all_species/pca_heterochromatin_erv_subspecies_pc1_pc2_maf_filter_20190708.pdf",
#   ggplot(
#   cbind(heterochromatin_erv_pca$x, indiv_df_heterochromatin_erv), 
#   aes(x=PC1, y=PC2, group=subspecies)) +
#   geom_point(aes(color=subspecies, shape=spectrum), size=5) +
#   scale_color_manual(values=subspecies_colors) +
#   xlab(generate_pc_axis_label(heterochromatin_erv_pca, 1)) +
#   ylab(generate_pc_axis_label(heterochromatin_erv_pca, 2)),
#   width = 7, height = 4)

# ggsave(
#   "./pca_all_species/pca_control_maternal_hotspots_chr16_subspecies_pc1_pc2_20190604.pdf",
#   ggplot(
#     cbind(control_maternal_hotspots_chr16_pca$x, indiv_df_control_maternal_hotspots),
#     aes(x=PC1, y=PC2, group=subspecies)) +
#     geom_point(aes(color=subspecies, shape=spectrum), size=5) +
#     scale_color_manual(values=subspecies_colors) +
#     xlab(generate_pc_axis_label(control_maternal_hotspots_chr16_pca, 1)) +
#     ylab(generate_pc_axis_label(control_maternal_hotspots_chr16_pca, 2)),
#   width = 7, height = 4)

# ggsave(
#   "./pca_all_species/pca_control_maternal_hotspots_chr16_subspecies_pc1_pc3_20190604.pdf",
#   ggplot(
#     cbind(control_maternal_hotspots_chr16_pca$x, indiv_df_control_maternal_hotspots),
#     aes(x=PC1, y=PC3, group=subspecies)) +
#     geom_point(aes(color=subspecies, shape=spectrum), size=5) +
#     scale_color_manual(values=subspecies_colors) +
#     xlab(generate_pc_axis_label(control_maternal_hotspots_pca, 1)) +
#     ylab(generate_pc_axis_label(control_maternal_hotspots_pca, 3)),
#   width = 7, height = 4)

# ggsave(
#   "./pca_all_species/pca_control_maternal_hotspots_subspecies_maf_filter_pc1_pc2_20190827.pdf",
#   ggplot(
#     cbind(control_maternal_hotspots_pca$x, indiv_df_control_maternal_hotspots),
#     aes(x=PC1, y=PC2, group=subspecies)) +
#     geom_point(aes(color=subspecies, shape=spectrum), size=5) +
#     scale_color_manual(values=subspecies_colors) +
#     xlab(generate_pc_axis_label(control_maternal_hotspots_pca, 1)) +
#     ylab(generate_pc_axis_label(control_maternal_hotspots_pca, 2)),
#   width = 7, height = 4)

# generate_pc_loading_heatmap(control_maternal_hotspots_pca$rotation[,2], "Control vs. Maternal Hotspots PC1")

# ggsave(
#   "./heatmaps/control_recom_pc1_20190603.pdf",
#   generate_pc_loading_heatmap(control_recom_pca$rotation[,1], "Control vs. Recombination PC1")
#   )

# ggsave(
#   "./heatmaps/control_recom_pc2_20190603.pdf",
#   generate_pc_loading_heatmap(control_recom_pca$rotation[,2], "Control vs. Recombination PC2")
# )

# ggsave(
#   "./heatmaps/control_recom_pc3_20190603.pdf",
#   generate_pc_loading_heatmap(control_recom_pca$rotation[,3], "Control vs. Recombination PC3")
# )

# ggsave(
#   "./heatmaps/control_recom_chimp_pc1_20190508.pdf",
#   generate_pc_loading_heatmap(control_recom_chimp_pca$rotation[,1], "Control vs. Chimp Recombination PC1")
# )

# ggsave(
#   "./heatmaps/control_recom_chimp_pc2_20190508.pdf",
#   generate_pc_loading_heatmap(control_recom_chimp_pca$rotation[,2], "Control vs. Chimp Recombination PC2")
# )


# ggsave(
#   "./heatmaps/control_recom_chimp_pc3_20190508.pdf",
#   generate_pc_loading_heatmap(control_recom_chimp_pca$rotation[,3], "Control vs. Chimp Recombination PC3")
# )

# ggsave(
#   "./heatmaps/control_early_late_rep_pca_pc1_20190603.pdf",
#   generate_pc_loading_heatmap(control_early_late_rep_pca$rotation[,1], "Control & Rep. Timing PC1")
# )

ggsave(
  "./heatmaps/control_early_late_rep_pca_pc2_20190827.pdf",
  generate_pc_loading_heatmap(control_early_late_rep_pca$rotation[,2], "Control & Rep. Timing PC2")
)

# pdf("./heatmaps/homo_indiv_pca_reweight_pc1.pdf")
# generate_pc_loading_heatmap(pca_by_species_repetitive_split_reweight$Homo$rotation[,1], "All compartments Homo reweight PC1")
# dev.off()

# pdf("./heatmaps/homo_indiv_pca_reweight_pc2.pdf")
# generate_pc_loading_heatmap(pca_by_species_repetitive_split_reweight$Homo$rotation[,2], "All compartments Homo reweight PC2")
# dev.off()

# pdf("./heatmaps/homo_indiv_pca_reweight_pc3.pdf")
# generate_pc_loading_heatmap(pca_by_species_repetitive_split_reweight$Homo$rotation[,3], "All compartments Homo reweight PC3")
# dev.off()

# pdf("./heatmaps/control_pca_pc1_20190603.pdf")
# generate_pc_loading_heatmap(control_pca$rotation[,1], "Control PCA PC1")
# dev.off()

# pdf("./heatmaps/control_pca_pc2_20190603.pdf")
# generate_pc_loading_heatmap(control_pca$rotation[,2], "Control PCA PC2")
# dev.off()

# A lot of this stuff doesn't require randomization. I haven't edited it yet.
early_rep_species_spectra <-
  load_species_3mer_spectrum("hg18_replication_timing_q3", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
late_rep_species_spectra <-
  load_species_3mer_spectrum("hg18_replication_timing_q0", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
early_rep_nonrepetitive_species_spectra <-
  load_species_3mer_spectrum("hg18_replication_timing_q3_nonrepetitive", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
late_rep_nonrepetitive_species_spectra <-
  load_species_3mer_spectrum("hg18_replication_timing_q0_nonrepetitive", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
early_rep_repetitive_species_spectra <-
  load_species_3mer_spectrum("hg18_replication_timing_q3_repetitive", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
late_rep_repetitive_species_spectra <-
  load_species_3mer_spectrum("hg18_replication_timing_q0_repetitive", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
erv_species_spectra <-
  load_species_3mer_spectrum("hg18_all_erv", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
control_species_spectra <-
  load_species_3mer_spectrum("hg18_control", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
# pericentromere_species_spectra <-
#   load_species_3mer_spectrum("hg18_pericentromeric_regions")
# subtelomere_species_spectra <-
#   load_species_3mer_spectrum("hg18_subtelomeric_regions")
# alu_species_spectra <-
#   load_species_3mer_spectrum("hg18_all_alu")
heterochromatin_species_spectra <-
  load_species_3mer_spectrum("hg18_nonrepetitive_chromHMM_heterochromatin", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
erv_hmc_low_species_spectra <-
  load_species_3mer_spectrum("hg18_all_erv_hmc_low", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
erv_hmc_high_species_spectra <-
  load_species_3mer_spectrum("hg18_all_erv_hmc_high", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
# maternal_hotspots_species_spectra <-
#   load_species_3mer_spectrum("hg18_maternal_hotspots_zero_indexed_mb", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
# old_control_species_spectra <- 
#   load_species_3mer_spectrum("hg18_nonrepetitive_neutral")

early_rep_species_spectra_reweight_control <-
  reweight_species_spectra(early_rep_species_spectra, 
                           early_rep_trinuc_composition,
                           control_trinuc_composition)
late_rep_species_spectra_reweight_control <-
  reweight_species_spectra(late_rep_species_spectra, 
                           late_rep_trinuc_composition,
                           control_trinuc_composition)
early_rep_nonrepetitive_species_spectra_reweight_control <-
  reweight_species_spectra(early_rep_nonrepetitive_species_spectra, 
                           early_rep_nonrepetitive_trinuc_composition,
                           control_trinuc_composition)
late_rep_nonrepetitive_species_spectra_reweight_control <-
  reweight_species_spectra(late_rep_nonrepetitive_species_spectra, 
                           late_rep_nonrepetitive_trinuc_composition,
                           control_trinuc_composition)
early_rep_repetitive_species_spectra_reweight_control <-
  reweight_species_spectra(early_rep_repetitive_species_spectra, 
                           early_rep_repetitive_trinuc_composition,
                           control_trinuc_composition)
late_rep_repetitive_species_spectra_reweight_control <-
  reweight_species_spectra(late_rep_repetitive_species_spectra, 
                           late_rep_repetitive_trinuc_composition,
                           control_trinuc_composition)
erv_species_spectra_reweight_control <-
  reweight_species_spectra(erv_species_spectra, 
                           erv_trinuc_composition,
                           control_trinuc_composition)
heterochromatin_species_spectra_reweight_control <-
  reweight_species_spectra(heterochromatin_species_spectra, 
                           heterochromatin_trinuc_composition,
                           control_trinuc_composition)
# pericentromere_species_spectra_reweight_control <-
#   reweight_species_spectra(pericentromere_species_spectra, 
#                            pericentromere_trinuc_composition,
#                            control_trinuc_composition)
# subtelomere_species_spectra_reweight_control <-
#   reweight_species_spectra(subtelomere_species_spectra, 
#                            subtelomere_trinuc_composition,
#                            control_trinuc_composition)
# alu_species_spectra_reweight_control <-
#   reweight_species_spectra(alu_species_spectra, 
#                            alu_trinuc_composition,
#                            control_trinuc_composition)
control_species_spectra_reweight_control <-
  reweight_species_spectra(control_species_spectra, 
                           control_trinuc_composition,
                           control_trinuc_composition)
erv_hmc_low_species_spectra_reweight_control <-
  reweight_species_spectra(erv_hmc_low_species_spectra, 
                           erv_hmc_low_trinuc_composition,
                           control_trinuc_composition)
erv_hmc_high_species_spectra_reweight_control <-
  reweight_species_spectra(erv_hmc_high_species_spectra, 
                           erv_hmc_high_trinuc_composition,
                           control_trinuc_composition)
# maternal_hotspots_species_spectra_reweight_control <-
#   reweight_species_spectra(maternal_hotspots_species_spectra, 
#                            maternal_hotspots_zero_indexed_trinuc_composition,
#                            control_trinuc_composition)

erv_hmc_fraction_for_ggplot <-
  data.frame(
    hmc_neg = unlist(erv_hmc_low_species_spectra_reweight_control),
    hmc_pos = unlist(erv_hmc_high_species_spectra_reweight_control),
    s = factor(rep(species, each=96)),
    mut = factor(collapsed_trinuc_mutations)
  )

ggsave(
  "./scatterplot_erv_hmc_mut_type_by_species_20190830.pdf",
  ggplot(erv_hmc_fraction_for_ggplot, aes(x=log(hmc_neg), y=log(hmc_pos), col=s, label=mut)) +
    geom_point() +
    geom_text(aes(label=ifelse(log(hmc_pos/hmc_neg)>0.4,as.character(mut),'')),hjust=0, vjust=0)
)

# Generating ERV HMC+ v HMC- p vals
erv_hmc_high_species_spectra_reweight_erv_hmc_low <-
  reweight_species_spectra(erv_hmc_high_species_spectra,
                           erv_hmc_high_trinuc_composition,
                           erv_hmc_low_trinuc_composition, for_chi_sq = T)
erv_hmc_low_species_spectra_reweight_erv_hmc_high <-
  reweight_species_spectra(erv_hmc_low_species_spectra,
                           erv_hmc_low_trinuc_composition,
                           erv_hmc_high_trinuc_composition, for_chi_sq = T)

write.csv(
  as.data.frame(
    sapply(
      species,
      function(s)
        chisq.test(
          matrix(
            c(
              sum(erv_hmc_high_species_spectra_reweight_erv_hmc_low[c("ACG.G", "CCG.G", "GCG.G", "TCG.G"), s]),
              sum(erv_hmc_low_species_spectra_reweight_erv_hmc_high[c("ACG.G", "CCG.G", "GCG.G", "TCG.G"), s]),
              sum(erv_hmc_high_species_spectra_reweight_erv_hmc_low[, s]) - 
                sum(erv_hmc_high_species_spectra_reweight_erv_hmc_low[c("ACG.G", "CCG.G", "GCG.G", "TCG.G"), s]),
              sum(erv_hmc_low_species_spectra_reweight_erv_hmc_high[, s]) -
                sum(erv_hmc_low_species_spectra_reweight_erv_hmc_high[c("ACG.G", "CCG.G", "GCG.G", "TCG.G"), s])
            ), nrow=2, ncol=2
          )
        )$p.value
    )
  ),
  "./cg_gg_p_vals_erv_hmc_high_to_low.csv"
)

ggsave(
  "./heatmaps/all_species_rep_timing_20190905.pdf",
  generate_heatmap_plot_multiple_species(late_rep_species_spectra_reweight_control,
                                         early_rep_species_spectra_reweight_control)
)

ggsave(
  "./heatmaps/all_species_erv_heterochromatin_20190913.pdf",
  generate_heatmap_plot_multiple_species(erv_species_spectra_reweight_control,
                                         heterochromatin_species_spectra_reweight_control)
)

ggsave(
  "./heatmaps/erv_hmc_high_vs_low_Homo.pdf",
  generate_heatmap_plot_single_species(
    erv_hmc_high_species_spectra_reweight_control,
    erv_hmc_low_species_spectra_reweight_control,
    "Homo"
  )
)

ggsave(
  "./heatmaps/erv_hmc_high_vs_low_Homo.pdf",
  generate_heatmap_plot_single_species(
    erv_hmc_high_species_spectra_reweight_control,
    erv_hmc_low_species_spectra_reweight_control,
    "Homo"
  )
)




ggsave(
  "./heatmaps/erv_hmc_high_vs_low_Pan_troglodytes.pdf",
  generate_heatmap_plot_single_species(
    erv_hmc_high_species_spectra_reweight_control,
    erv_hmc_low_species_spectra_reweight_control,
    "Pan_troglodytes"
  )
)
ggsave(
  "./heatmaps/erv_hmc_high_vs_low_Pan_paniscus.pdf",
  generate_heatmap_plot_single_species(
    erv_hmc_high_species_spectra_reweight_control,
    erv_hmc_low_species_spectra_reweight_control,
    "Pan_paniscus"
  )
)
ggsave(
  "./heatmaps/heterochromatin_vs_erv_Pongo_abelii.pdf",
  generate_heatmap_plot_single_species(
    erv_species_spectra_reweight_control,
    heterochromatin_species_spectra_reweight_control,
    "Pongo_abelii"
  )
)
ggsave(
  "./heatmaps/heterochromatin_vs_erv_Homo.pdf",
  generate_heatmap_plot_single_species(
    erv_species_spectra_reweight_control,
    heterochromatin_species_spectra_reweight_control,
    "Homo"
  )
)
ggsave(
  "./heatmaps/heterochromatin_vs_erv_Gorilla.pdf",
  generate_heatmap_plot_single_species(
    erv_species_spectra_reweight_control,
    heterochromatin_species_spectra_reweight_control,
    "Gorilla"
  )
)

ggsave(
  "./heatmaps/erv_hmc_high_vs_low_Pongo_abelii.pdf",
  generate_heatmap_plot_single_species(
    erv_hmc_high_species_spectra_reweight_control,
    erv_hmc_low_species_spectra_reweight_control,
    "Pongo_abelii"
  )
)
ggsave(
  "./heatmaps/erv_hmc_high_vs_low_Pongo_pygmaeus.pdf",
  generate_heatmap_plot_single_species(
    erv_hmc_high_species_spectra_reweight_control,
    erv_hmc_low_species_spectra_reweight_control,
    "Pongo_pygmaeus"
  )
)
ggsave(
  "./heatmaps/maternal_hotspots_vs_control_Homo.pdf",
  generate_heatmap_plot_single_species(
    maternal_hotspots_species_spectra_reweight_control,
    control_species_spectra_reweight_control,
    "Homo"
  )
)
ggsave(
  "./heatmaps/late_vs_early_rep_Homo.pdf",
  generate_heatmap_plot_single_species(
    late_rep_species_spectra_reweight_control,
    early_rep_species_spectra_reweight_control,
    "Homo"
  )
)
ggsave(
  "./heatmaps/late_vs_early_rep_Gorilla.pdf",
  generate_heatmap_plot_single_species(
    late_rep_species_spectra_reweight_control,
    early_rep_species_spectra_reweight_control,
    "Gorilla"
  )
)
ggsave(
  "./heatmaps/late_vs_early_rep_Pan_troglodytes.pdf",
  generate_heatmap_plot_single_species(
    late_rep_species_spectra_reweight_control,
    early_rep_species_spectra_reweight_control,
    "Pan_troglodytes"
  )
)
ggsave(
  "./heatmaps/late_vs_early_rep_Pongo_abelii.pdf",
  generate_heatmap_plot_single_species(
    late_rep_species_spectra_reweight_control,
    early_rep_species_spectra_reweight_control,
    "Pongo_abelii"
  )
)

ggsave(
  "./heatmaps/late_vs_early_rep_repetitive_Gorilla.pdf",
  generate_heatmap_plot_single_species(
    late_rep_repetitive_species_spectra_reweight_control,
    early_rep_repetitive_species_spectra_reweight_control,
    "Gorilla"
  )
)

ggsave(
  "./heatmaps/late_vs_early_rep_nonrepetitive_Gorilla.pdf",
  generate_heatmap_plot_single_species(
    late_rep_nonrepetitive_species_spectra_reweight_control,
    early_rep_nonrepetitive_species_spectra_reweight_control,
    "Gorilla"
  )
)



ggsave(
  "./heatmaps/maternal_hotspots_vs_control_Pan_troglodytes.pdf",
  generate_heatmap_plot_single_species(
    maternal_hotspots_species_spectra_reweight_control,
    control_species_spectra_reweight_control,
    "Pan_troglodytes"
  )
)
ggsave(
  "./heatmaps/maternal_hotspots_vs_control_Pongo_abelii.pdf",
  generate_heatmap_plot_single_species(
    maternal_hotspots_species_spectra_reweight_control,
    control_species_spectra_reweight_control,
    "Pongo_abelii"
  )
)

ggsave(
  "./heatmaps/Homo_erv_vs_control_20190627.pdf",
  generate_heatmap_plot_single_species(
    erv_species_spectra_reweight_control,
    control_species_spectra_reweight_control,
    "Pongo_abelii"
  )
)

ggsave(
  "./heatmaps/Homo_late_rep_nr_vs_early_rep_nr_20190627.pdf",
  generate_heatmap_plot_single_species(
    late_rep_nonrepetitive_species_spectra_reweight_control,
    early_rep_nonrepetitive_species_spectra_reweight_control,
    "Homo"
  )
)

ggsave(
  "./heatmaps/Pan_troglodytes_late_rep_nr_vs_early_rep_nr_20190627.pdf",
  generate_heatmap_plot_single_species(
    late_rep_nonrepetitive_species_spectra_reweight_control,
    early_rep_nonrepetitive_species_spectra_reweight_control,
    "Pan_troglodytes"
  )
)
ggsave(
  "./heatmaps/Gorilla_late_rep_nr_vs_early_rep_nr_20190627.pdf",
  generate_heatmap_plot_single_species(
    late_rep_nonrepetitive_species_spectra_reweight_control,
    early_rep_nonrepetitive_species_spectra_reweight_control,
    "Gorilla"
  )
)

generate_heatmap_plot_single_species(
  maternal_hotspots_species_spectra_reweight_control,
  control_species_spectra_reweight_control,
  "Pan_troglodytes"
)

generate_heatmap_plot_single_species(
  maternal_hotspots_species_spectra_reweight_control,
  control_species_spectra_reweight_control,
  "Pongo_abelii"
)

log_odds_spectra_correlation_species(
  erv_species_spectra_reweight_control,
  control_species_spectra_reweight_control
)

# ggsave(
#   "./correlation_plots/lr_erv_vs_control_corr_spp_20190627.pdf",
#   log_odds_spectra_correlation_species_heatmap(
#     erv_species_spectra_reweight_control,
#     control_species_spectra_reweight_control
#   )
# )
# ggsave(
#   "./correlation_plots/corr_erv_vs_control_lr_spp_20190627.pdf",
#   log_odds_species_correlation_spectra_heatmap(
#     erv_species_spectra_reweight_control,
#     control_species_spectra_reweight_control
#   )
# )
# ggsave(
#   "./correlation_plots/lr_late_rep_nr_vs_early_rep_nr_corr_spp_20190627.pdf",
#   log_odds_spectra_correlation_species_heatmap(
#     late_rep_nonrepetitive_species_spectra_reweight_control,
#     early_rep_nonrepetitive_species_spectra_reweight_control
#   )
# )
# ggsave(
#   "./correlation_plots/corr_late_rep_nr_vs_early_rep_nr_lr_spp_20190627.pdf",
#   log_odds_species_correlation_spectra_heatmap(
#     late_rep_nonrepetitive_species_spectra_reweight_control,
#     early_rep_nonrepetitive_species_spectra_reweight_control
#   )
# )

write.csv(
  log_odds_spectra_correlation_species(
    late_rep_species_spectra_reweight_control,
    early_rep_species_spectra_reweight_control
  ),
  "./corr_values_log_odds_rep_timing_corr_species.csv"
)

# log_odds_spectra_correlation_species(
#   late_rep_species_spectra_reweight_control,
#   early_rep_species_spectra_reweight_control
# )


# pdf("./heatmaps/alu_vs_control_reweight_control.pdf")
# generate_heatmap_plot_single_species(
#   alu_species_spectra_reweight_control,
#   control_species_spectra_reweight_control,
#   "Homo"
# )
# dev.off()

# pdf("./heatmaps/subtelomere_vs_control_reweight_control.pdf")
# generate_heatmap_plot_single_species(
#   subtelomere_species_spectra_reweight_control,
#   control_species_spectra_reweight_control,
#   "Homo"
# )
# dev.off()

# pdf("./heatmaps/pericentromere_vs_control_reweight_control.pdf")
# generate_heatmap_plot_single_species(
#   pericentromere_species_spectra_reweight_control,
#   control_species_spectra_reweight_control,
#   "Homo"
# )
# dev.off()

# generate_heatmap_plot_single_species(
#   pericentromere_species_spectra_reweight_control,
#   late_rep_repetitive_species_spectra_reweight_control,
#   "Homo"
# )

# generate_heatmap_plot_single_species(
#   subtelomere_species_spectra_reweight_control,
#   early_rep_repetitive_species_spectra_reweight_control,
#   "Homo"
# )

# pdf("./heatmaps/early_rep_nonrepetitive_vs_late_rep_nonrepetitive_reweight_control.pdf")
# generate_heatmap_plot_single_species(
#   early_rep_nonrepetitive_species_spectra_reweight_control,
#   late_rep_nonrepetitive_species_spectra_reweight_control,
#   "Homo"
# )
# dev.off()

ggsave(
  "./heatmaps/erv_vs_control_chisq_homo_20190619.pdf",
  generate_heatmap_plot_single_species_chisq(
    erv_species_spectra,
    erv_trinuc_composition,
    control_species_spectra,
    control_trinuc_composition,
    "Homo"
  )
)



ggsave(
  "./heatmaps/erv_vs_control_species_20190508.pdf",
  generate_heatmap_plot_single_species(
    erv_species_spectra_reweight_control,
    control_species_spectra_reweight_control,
    "Homo"
  )
  )

ggsave(
  "./heatmaps/erv_vs_heterochromatin_species_20190508.pdf",
  generate_heatmap_plot_single_species(
    erv_species_spectra_reweight_control,
    heterochromatin_species_spectra_reweight_control,
    "Homo"
  )
)

# pdf("./heatmaps/early_vs_late_replication_timing_nonrepetitive_nonnormalized.pdf")
generate_pc_loading_heatmap(
  log(early_rep_nonrepetitive_species_spectra_fraction_homo / late_rep_nonrepetitive_species_spectra_fraction_homo), 
  "Non-normaized early vs. late replication nonrepetitive")
# dev.off()

generate_pc_loading_heatmap(
  log(early_rep_nonrepetitive_species_spectra_reweight$Homo / late_rep_nonrepetitive_species_spectra$Homo), 
  "Non-normaized early vs. late replication nonrepetitive")


generate_pc_loading_heatmap(
  log((early_rep_nonrepetitive_species_spectra_reweight$Homo /
         sum(early_rep_nonrepetitive_species_spectra_reweight$Homo))
      / (late_rep_nonrepetitive_species_spectra$Homo /
           sum(late_rep_nonrepetitive_species_spectra$Homo))), 
  "Reweighted early vs. late replication nonrepetitive")


compare_early_late_trinuc <-
  data.frame(
    logodds = unlist(log(early_rep_nonrepetitive_trinuc_composition / late_rep_nonrepetitive_trinuc_composition)),
    five_prime_and_center = as.factor(substr(names(early_rep_nonrepetitive_trinuc_composition), 1, 2)),
    three_prime = as.factor(substr(names(early_rep_nonrepetitive_trinuc_composition), 3, 3)))

pdf("./heatmaps/early_vs_late_replication_timing_nonrepetitive_trinuc_content.pdf")
ggplot(compare_early_late_trinuc, aes(three_prime, five_prime_and_center))+
  geom_tile(aes(fill=logodds), color="white") +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = 0)
dev.off()

pdf("./heatmaps/early_rep_nonrepetitive_vs_late_rep_nonrepetitive.pdf")
generate_heatmap_plot_single_species(
  early_rep_nonrepetitive_species_spectra_norm,
  late_rep_nonrepetitive_species_spectra_norm,
  s="Homo"
)
dev.off()

# pdf("./heatmaps/early_rep_nonrepetitive_vs_late_rep_nonrepetitive.pdf")
generate_heatmap_plot_single_species( # ugh this isn't working
  early_rep_nonrepetitive_species_spectra_reweight,
  late_rep_nonrepetitive_species_spectra,
  s="Homo"
)
# dev.off()

pdf("./heatmaps/early_rep_repetitive_vs_late_rep_repetitive.pdf")
generate_heatmap_plot_single_species(
  early_rep_repetitive_species_spectra_norm,
  late_rep_repetitive_species_spectra_norm,
  s="Homo"
) 
dev.off()

# Looking at what's up with Donald

pan_t_v_signature <-
  sapply(
    names(control_spectra_reweight),
    function(x)
      mean(
        subset(
          control_spectra_reweight[, x],
          startsWith(rownames(control_spectra_reweight), "Pan_troglodytes_verus") &
            !endsWith(rownames(control_spectra_reweight), "Donald")
        )
      )
  )
pan_t_t_signature <-
  sapply(
    names(control_spectra_reweight),
    function(x)
      mean(
        subset(
          control_spectra_reweight[, x],
          startsWith(rownames(control_spectra_reweight), "Pan_troglodytes_troglodytes")
        )
      )
  )
dist(
  rbind(
    pan_t_v_signature, 
    pan_t_t_signature
    )
  )
dist(
  rbind(
    pan_t_v_signature, 
    subset(control_spectra_reweight, 
           endsWith(rownames(control_spectra_reweight), "Donald")
           )
    )
  )
dist(
  rbind(
    pan_t_t_signature, 
    subset(control_spectra_reweight, 
           endsWith(rownames(control_spectra_reweight), "Donald")
    )
  )
)

# ggsave(
#   "./violin_plot_control_20190526.pdf",
#   plot_indiv_compartment_chisq(
#     control_spectra,
#     control_trinuc_composition, 
#     control_spectra,
#     control_trinuc_composition,
#     "Compare control spectra"
#   )
# )


# I AM HERE!

# library(umap)
# erv_and_control_umap <- umap(erv_and_control_spectra_norm)
# colnames(erv_and_control_umap$layout) <- c("x", "y")
# plot(erv_and_control_umap$layout, col=indiv_df_erv_and_control_spectra$col)

# pdf("./erv_vs_control_umap_20190109.pdf", width=6, height=4.5)
# ggplot(
#   cbind(erv_and_control_umap$layout, indiv_df_erv_and_control_spectra), 
#   aes(x=x, y=y, group=subspecies)) +
#   geom_point(aes(color=subspecies, shape=spectrum), size=2.5)
# dev.off()
