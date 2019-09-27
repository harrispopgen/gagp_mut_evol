# Analyze mutational signatures per individual (lots of PCAs)

# Set working directory here
# setwd("~/Documents/Harris_project_directory")

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

ncnr_trinuc_composition <- read.table("./hg18_ncnr_3mer_content.txt", header = T)
erv_trinuc_composition <- read.table("./hg18_all_erv_3mer_content.txt", header = T)
heterochromatin_trinuc_composition <- read.table("./hg18_nonrepetitive_chromHMM_heterochromatin_3mer_content.txt", header = T)
line_trinuc_composition <- read.table("./hg18_all_line_3mer_content.txt", header = T)
early_rep_nonrepetitive_trinuc_composition <- read.table("./hg18_replication_timing_q3_nonrepetitive_3mer_content.txt", header = T)
early_rep_repetitive_trinuc_composition <- read.table("./hg18_replication_timing_q3_repetitive_3mer_content.txt", header = T)
late_rep_nonrepetitive_trinuc_composition <- read.table("./hg18_replication_timing_q0_nonrepetitive_3mer_content.txt", header = T)
late_rep_repetitive_trinuc_composition <- read.table("./hg18_replication_timing_q0_repetitive_3mer_content.txt", header = T)

ncnr_spectra <- load_indiv_3mer_spectrum("hg18_ncnr", nmer_dir = "./nmer_mutation_counts/")
erv_spectra <- load_indiv_3mer_spectrum("hg18_all_erv", nmer_dir = "./nmer_mutation_counts/")
heterochromatin_spectra <- load_indiv_3mer_spectrum("hg18_nonrepetitive_chromHMM_heterochromatin", nmer_dir = "./nmer_mutation_counts/")
line_spectra <- load_indiv_3mer_spectrum("hg18_all_line", nmer_dir = "./nmer_mutation_counts/")
early_rep_repetitive_spectra <- load_indiv_3mer_spectrum("hg18_replication_timing_q3_repetitive", nmer_dir = "./nmer_mutation_counts/")
early_rep_nonrepetitive_spectra <- load_indiv_3mer_spectrum("hg18_replication_timing_q3_nonrepetitive", nmer_dir = "./nmer_mutation_counts/")
late_rep_repetitive_spectra <- load_indiv_3mer_spectrum("hg18_replication_timing_q0_repetitive", nmer_dir = "./nmer_mutation_counts/")
late_rep_nonrepetitive_spectra <- load_indiv_3mer_spectrum("hg18_replication_timing_q0_nonrepetitive", nmer_dir = "./nmer_mutation_counts/")

ncnr_spectra_reweight <- reweight_indiv_spectra(ncnr_spectra, ncnr_trinuc_composition, ncnr_trinuc_composition)
erv_spectra_reweight <- reweight_indiv_spectra(erv_spectra, erv_trinuc_composition, ncnr_trinuc_composition)
heterochromatin_spectra_reweight <- reweight_indiv_spectra(heterochromatin_spectra, heterochromatin_trinuc_composition, ncnr_trinuc_composition)
line_spectra_reweight <- reweight_indiv_spectra(line_spectra, line_trinuc_composition, ncnr_trinuc_composition)
early_rep_repetitive_spectra_reweight <- reweight_indiv_spectra(early_rep_repetitive_spectra, early_rep_repetitive_trinuc_composition, ncnr_trinuc_composition)
early_rep_nonrepetitive_spectra_reweight <- reweight_indiv_spectra(early_rep_nonrepetitive_spectra, early_rep_nonrepetitive_trinuc_composition, ncnr_trinuc_composition)
late_rep_repetitive_spectra_reweight <- reweight_indiv_spectra(late_rep_repetitive_spectra, late_rep_repetitive_trinuc_composition, ncnr_trinuc_composition)
late_rep_nonrepetitive_spectra_reweight <- reweight_indiv_spectra(late_rep_nonrepetitive_spectra, late_rep_nonrepetitive_trinuc_composition, ncnr_trinuc_composition)

spectra_by_species_reweight <-
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
names(spectra_by_species_reweight) <- species

indiv_df_by_species <-
  lapply(
    species,
    function(s){
      n_indiv = nrow(subset(control_spectra, startsWith(rownames(control_spectra), s)))
      return(
        data.frame(
          "indiv" = rownames(spectra_by_species_reweight[[s]]),
          "subspecies" = sub("-.*", "", rownames(spectra_by_species_reweight[[s]])),
          "spectrum" = c(rep("ncnr", n_indiv),
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
names(indiv_df_by_species) <- species

pca_by_species <-
  lapply(
    spectra_by_species_reweight,
    function(spectra)
      prcomp(
        spectra,
        scale. = T
      )
  )
names(pca_by_species) <- species

ggsave(
  "./images/all_spectrum_pca_Homo_pc1_pc2.pdf",
  ggplot(
    cbind(pca_by_species$Homo$x, indiv_df_by_species$Homo),
    aes(x=PC1, y=PC2, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[4]) +
    scale_fill_manual(values=subspecies_colors[4]) +
    xlab(generate_pc_axis_label(pca_by_species$Homo, 1)) +
    ylab(generate_pc_axis_label(pca_by_species$Homo, 2)),
  width = 7, height = 4)

ggsave(
  "./images/all_spectrum_pca_Homo_pc1_pc3.pdf",
  ggplot(
    cbind(pca_by_species$Homo$x, indiv_df_by_species$Homo),
    aes(x=PC1, y=PC3, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[4]) +
    scale_fill_manual(values=subspecies_colors[4]) +
    xlab(generate_pc_axis_label(pca_by_species$Homo, 1)) +
    ylab(generate_pc_axis_label(pca_by_species$Homo, 3)),
  width = 7, height = 4)

ggsave(
  "./images/all_spectrum_pca_Pan_troglodytes_pc1_pc2.pdf",
  ggplot(
    cbind(pca_by_species$Pan_troglodytes$x, indiv_df_by_species$Pan_troglodytes),
    aes(x=PC1, y=PC2, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[c(6:9)]) +
    scale_fill_manual(values=subspecies_colors[c(6:9)]) +
    xlab(generate_pc_axis_label(pca_by_species$Pan_troglodytes, 1)) +
    ylab(generate_pc_axis_label(pca_by_species$Pan_troglodytes, 2)),
  width = 7, height = 4)

ggsave(
  "./images/all_spectrum_pca_Pan_troglodytes_pc1_pc3.pdf",
  ggplot(
    cbind(pca_by_species$Pan_troglodytes$x, indiv_df_by_species$Pan_troglodytes),
    aes(x=PC1, y=PC3, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[c(6:9)]) +
    scale_fill_manual(values=subspecies_colors[c(6:9)]) +
    xlab(generate_pc_axis_label(pca_by_species$Pan_troglodytes, 1)) +
    ylab(generate_pc_axis_label(pca_by_species$Pan_troglodytes, 3)),
  width = 7, height = 4)

ggsave(
  "./images/all_spectrum_pca_Pongo_abelii_pc1_pc2.pdf",
  ggplot(
    cbind(pca_by_species$Pongo_abelii$x, indiv_df_by_species$Pongo_abelii),
    aes(x=PC1, y=PC2, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[c(10)]) +
    scale_fill_manual(values=subspecies_colors[c(10)]) +
    xlab(generate_pc_axis_label(pca_by_species$Pongo_abelii, 1)) +
    ylab(generate_pc_axis_label(pca_by_species$Pongo_abelii, 2)),
  width = 7, height = 4)

ggsave(
  "./images/all_spectrum_pca_Pongo_abelii_pc1_pc3.pdf",
  ggplot(
    cbind(pca_by_species$Pongo_abelii$x, indiv_df_by_species$Pongo_abelii),
    aes(x=PC1, y=PC3, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(16, 17, 2, 3, 18, 15, 0, 4))+
    scale_color_manual(values=subspecies_colors[c(10)]) +
    scale_fill_manual(values=subspecies_colors[c(10)]) +
    xlab(generate_pc_axis_label(pca_by_species$Pongo_abelii, 1)) +
    ylab(generate_pc_axis_label(pca_by_species$Pongo_abelii, 3)),
  width = 7, height = 4)