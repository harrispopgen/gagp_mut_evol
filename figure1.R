# Analyze mutational signatures per individual (lots of PCAs)

# Set working directory here
# setwd("~/Documents/Harris_project_directory")

# Need to comment which figure panels

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
early_rep_trinuc_composition <- read.table("./hg18_replication_timing_q3_3mer_content.txt", header = T)
late_rep_trinuc_composition <- read.table("./hg18_replication_timing_q0_3mer_content.txt", header = T)

ncnr_spectra <- load_indiv_3mer_spectrum("hg18_ncnr", nmer_dir = "./nmer_mutation_counts/")
early_rep_spectra <- load_indiv_3mer_spectrum("hg18_replication_timing_q3", nmer_dir = "./nmer_mutation_counts/")
late_rep_spectra <- load_indiv_3mer_spectrum("hg18_replication_timing_q0", nmer_dir = "./nmer_mutation_counts/")

ncnr_spectra_reweight <- reweight_indiv_spectra(ncnr_spectra, ncnr_trinuc_composition, ncnr_trinuc_composition)
early_rep_spectra_reweight <- reweight_indiv_spectra(early_rep_spectra, early_rep_trinuc_composition, ncnr_trinuc_composition)
late_rep_spectra_reweight <- reweight_indiv_spectra(late_rep_spectra, late_rep_trinuc_composition, ncnr_trinuc_composition)

indiv_df <-
  data.frame(
    "indiv" = rownames(ncnr_spectra_reweight),
    "subspecies" = sub("-.*", "", rownames(ncnr_spectra_reweight))
  )
indiv_df$species <-
  sapply(indiv_df$subspecies,
         function(x)
           if(startsWith(x, "Gorilla")) "Gorilla" else 
             if(startsWith(x, "Pan_troglodytes")) "Pan_troglodytes" else x)
ncnr_pca <-
  prcomp(
    ncnr_spectra_reweight,
    scale. = T
  )

indiv_df_ncnr_rep_time <- 
  data.frame(
    "indiv" = 
      c(rownames(ncnr_spectra_reweight),
        rownames(early_rep_spectra_reweight),
        rownames(late_rep_spectra_reweight)),
    "subspecies"=
      c(
        sub("-.*", "", rownames(ncnr_spectra_reweight)),
        sub("-.*", "", rownames(early_rep_spectra_reweight)),
        sub("-.*", "", rownames(late_rep_spectra_reweight))
      ),
    "spectrum" = c(rep("ncnr", nrow(ncnr_spectra_reweight)),
                   rep("early_rep", nrow(early_rep_spectra_reweight)),
                   rep("late_rep", nrow(late_rep_spectra_reweight)))
  )
indiv_df_ncnr_rep_time$species <-
  sapply(indiv_df_ncnr_rep_time$subspecies,
         function(x)
           if(startsWith(x, "Gorilla")) "Gorilla" else 
             if(startsWith(x, "Pan_troglodytes")) "Pan_troglodytes" else x)

ncnr_early_late_rep_pca <-
  prcomp(
    rbind(
      ncnr_spectra_reweight,
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
        subset(ncnr_spectra_reweight, startsWith(rownames(ncnr_spectra_reweight), s)),
        subset(early_rep_spectra_reweight, startsWith(rownames(ncnr_spectra_reweight), s)),
        subset(late_rep_spectra_reweight, startsWith(rownames(ncnr_spectra_reweight), s))
      )
  )
names(spectra_by_species_replication_timing) <- species

ncnr_early_late_rep_pca_by_species <-
  lapply(
    species,
    function(s)
      prcomp(
        spectra_by_species_replication_timing[[s]],
        scale. = T
      )
  )
names(ncnr_early_late_rep_pca_by_species) <- species

indiv_df_by_species_ncnr_early_late_rep <-
  lapply(
    species,
    function(s){
      n_indiv = nrow(subset(ncnr_spectra, startsWith(rownames(ncnr_spectra), s)))
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
names(indiv_df_by_species_ncnr_early_late_rep) <- species

ggsave(
  "./images/rep_timing_pca_Homo_pc1_pc2.pdf",
  ggplot(
    cbind(ncnr_early_late_rep_pca_by_species$Homo$x, indiv_df_by_species_ncnr_early_late_rep$Homo),
    aes(x=PC1, y=PC2, group=subspecies)) +
    geom_point(aes(shape=spectrum, color=subspecies, fill = subspecies), size=5) +
    scale_shape_manual(values = c(17, 15, 16))+
    scale_color_manual(values=subspecies_colors[4]) +
    scale_fill_manual(values=subspecies_colors[4]) +
    xlab(generate_pc_axis_label(ncnr_early_late_rep_pca_by_species$Homo, 1)) +
    ylab(generate_pc_axis_label(ncnr_early_late_rep_pca_by_species$Homo, 2)),
  width = 7, height = 4)

ggsave("./images/pca_ncnr_pc1_pc2.pdf",
       ggplot(
         cbind(ncnr_pca$x, indiv_df),
         aes(x=PC1, y=PC2, group=subspecies)) +
         geom_point(aes(fill=subspecies), colour="black", pch=21, size=5) +
         scale_fill_manual(values=subspecies_colors) +
         xlab(generate_pc_axis_label(ncnr_pca, 1)) +
         ylab(generate_pc_axis_label(ncnr_pca, 2)),
       width=7, height=4)

ggsave(
  "./images/pca_ncnr_rep_timing_pc1_pc2.pdf",
  ggplot(
    cbind(ncnr_early_late_rep_pca$x, indiv_df_ncnr_rep_time), 
    aes(x=PC1, y=PC2, group=subspecies)) +
    geom_point(aes(color=subspecies, shape=spectrum), size=5) +
    scale_color_manual(values=subspecies_colors) +
    xlab(generate_pc_axis_label(ncnr_early_late_rep_pca, 1)) +
    ylab(generate_pc_axis_label(ncnr_early_late_rep_pca, 2)),
  width = 7, height = 4)

ggsave("./images/ridgeplot_subset_ncnr.pdf",
       plot_indiv_compartment_chisq(ncnr_spectra,
                                    ncnr_trinuc_composition,
                                    ncnr_spectra,
                                    ncnr_trinuc_composition,
                                    "NCNR Chi2 Distribution",
                                    c("Homo", "Pan_troglodytes", "Gorilla"))
)

early_rep_species_spectra <-
  load_species_3mer_spectrum("hg18_replication_timing_q3", nmer_dir = "./nmer_mutation_counts/")
late_rep_species_spectra <-
  load_species_3mer_spectrum("hg18_replication_timing_q0", nmer_dir = "./nmer_mutation_counts/")

early_rep_species_spectra_reweight_ncnr <-
  reweight_species_spectra(early_rep_species_spectra, 
                           early_rep_trinuc_composition,
                           ncnr_trinuc_composition)
late_rep_species_spectra_reweight_ncnr <-
  reweight_species_spectra(late_rep_species_spectra, 
                           late_rep_trinuc_composition,
                           ncnr_trinuc_composition)

ggsave(
  "./images/late_v_early_rep_heatmap_Homo.pdf",
  generate_heatmap_plot_single_species(
    late_rep_species_spectra_reweight_ncnr,
    early_rep_species_spectra_reweight_ncnr,
    "Homo"
  )
)

ggsave(
  "./images/all_species_late_v_early_rep_heatmap.pdf",
  generate_heatmap_plot_multiple_species(late_rep_species_spectra_reweight_ncnr,
                                         early_rep_species_spectra_reweight_ncnr)
)