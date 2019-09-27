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
erv_hmc_high_trinuc_composition <- read.table("./hg18_all_erv_hmc_high_3mer_content.txt", header = T)
erv_hmc_low_trinuc_composition <- read.table("./hg18_all_erv_hmc_low_3mer_content.txt", header = T)

# A lot of this stuff doesn't require randomization. I haven't edited it yet.
erv_species_spectra <-
  load_species_3mer_spectrum("hg18_all_erv", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
control_species_spectra <-
  load_species_3mer_spectrum("hg18_control", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
heterochromatin_species_spectra <-
  load_species_3mer_spectrum("hg18_nonrepetitive_chromHMM_heterochromatin", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
erv_hmc_low_species_spectra <-
  load_species_3mer_spectrum("hg18_all_erv_hmc_low", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")
erv_hmc_high_species_spectra <-
  load_species_3mer_spectrum("hg18_all_erv_hmc_high", nmer_dir = "./snp_data/randomized_nmer_mutation_counts_maf_filter_exclude_recurrent/")

erv_species_spectra_reweight_control <-
  reweight_species_spectra(erv_species_spectra, 
                           erv_trinuc_composition,
                           control_trinuc_composition)
heterochromatin_species_spectra_reweight_control <-
  reweight_species_spectra(heterochromatin_species_spectra, 
                           heterochromatin_trinuc_composition,
                           control_trinuc_composition)
erv_hmc_low_species_spectra_reweight_control <-
  reweight_species_spectra(erv_hmc_low_species_spectra, 
                           erv_hmc_low_trinuc_composition,
                           control_trinuc_composition)
erv_hmc_high_species_spectra_reweight_control <-
  reweight_species_spectra(erv_hmc_high_species_spectra, 
                           erv_hmc_high_trinuc_composition,
                           control_trinuc_composition)

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
