# Some code that I use pretty frequently for GAGP mutational signature analyses

library(ggplot2)
library(ggfortify)
library(ggridges)
library(scales)
library(grid)
library(gridExtra)
library(cowplot)
library(plyr)

setwd("~/Documents/Harris_project_directory")

chrs <- paste0("chr", c(1:22))
species <- c("Homo", "Pan_troglodytes", "Pan_paniscus", "Gorilla", "Pongo_abelii", "Pongo_pygmaeus")

collapsed_trinuc_mutations <- 
  paste0(
    rep(c("A", "C", "G", "T"), each=24),
    rep(c("A", "C"), each = 12, times=4),
    rep(c("A", "C", "G", "T"), each=3, times=8),
    rep(".", times=96),
    rep(
      c(rep(c("C", "G", "T"), times=4),
        rep(c("A", "G", "T"), times=4)),
      times=4
    )
  )

collapsed_trinuc_mutations_for_heatmap <-
  paste0(
    rep(c("A", "C", "G", "T"), each = 4, times = 6),
    rep(c("A", "C"), each = 48),
    rep(c("A", "C", "G", "T"), times = 24),
    rep(".", times = 96),
    rep(c("T", "C", "G", "T", "G", "A"), each = 16)
  )

five_prime_and_mut <- 
  factor(
    paste0(
      substr(collapsed_trinuc_mutations, 1, 2),
      substr(collapsed_trinuc_mutations, 4, 5)),
    levels = 
      unique(
        paste0(
          substr(collapsed_trinuc_mutations_for_heatmap, 1, 2),
          substr(collapsed_trinuc_mutations_for_heatmap, 4, 5)
        )
      )
  )

three_prime <- 
  factor(
    substr(collapsed_trinuc_mutations, 3, 3),
    levels = 
      unique(
        substr(collapsed_trinuc_mutations_for_heatmap, 3, 3)
      )
  )

generate_pc_loading_heatmap <- function(pca_loading, plot_title) {
  db = 
    data.frame(
      pca_loading = pca_loading,
      five_prime_and_mut = five_prime_and_mut,
      three_prime = three_prime
    )
  extreme_val = max(abs(min(pca_loading)), max(pca_loading))
  plot_list =
    lapply(
      rev(unique(substr(levels(db[, c("five_prime_and_mut")]), 2, 4))),
      function(mut)
        ggplot(subset(db, substr(db[, c("five_prime_and_mut")], 2, 4) == mut), aes(three_prime, five_prime_and_mut)) +
        geom_tile(aes(fill=subset(db[, c("pca_loading")], substr(db[, c("five_prime_and_mut")], 2, 4) == mut)), color = "white") +
        scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = 0, limits = c(-1*extreme_val, extreme_val)) +
        theme_minimal() +
        coord_equal() +
        ylab(sub("\\.", ">", mut)) +
        scale_y_discrete(labels=c("A", "C", "G", "T")) +
        theme(axis.text.x = element_blank(), 
              axis.ticks = element_blank(), 
              axis.title.x = element_blank(), 
              axis.title.y = element_text(angle = 0, vjust = 0.5),
              legend.position = "none",
              panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "white"),
              plot.margin =unit(c(0,0,0,0), "cm"))
    )
  plot_list[[1]] =
    plot_list[[1]] +
    ggtitle(plot_title) +
    theme(plot.title=element_text())
  plot_list[[6]] = 
    plot_list[[6]] + 
    theme(axis.text.x=element_text())
  for_legend = 
    ggplot(db, aes(three_prime, five_prime_and_mut)) +
    geom_tile(aes(fill=db[, c("pca_loading")]), color = "white") +
    scale_fill_gradient2(low = muted("blue"), 
                         mid = "white", 
                         high = muted("red"), 
                         midpoint = 0, 
                         limits = c(-1*extreme_val, extreme_val)) +
    theme(legend.title = element_blank())
  plot_list[[7]] <- cowplot::get_legend(for_legend)
  hlay = rbind(c(1,7),
                c(2,7),
                c(3,7),
                c(4,7),
                c(5,7),
                c(6,7))
  grid.arrange(grobs=plot_list, 
               ncol = 2, 
               layout_matrix = hlay,
               heights = unit(c(2.6, rep(2, 4), 2.25), "cm"))
}

load_indiv_3mer_spectrum = 
  function(compartment_name, nmer_dir = "~/Documents/Harris_project_directory/snp_data/nmer_mutation_counts/"){
    spectra_list = 
      lapply(species, function(s){
        filename = paste0(nmer_dir, s, "_", compartment_name, "_indiv_3mer_counts.txt")
        return(as.data.frame(read.table(filename, header = T, row.names = 1)))
      })
    names(spectra_list) <- species
    spectra = 
      rbind(
        spectra_list$Homo,
        spectra_list$Pan_troglodytes,
        spectra_list$Pan_paniscus,
        spectra_list$Gorilla,
        spectra_list$Pongo_abelii,
        spectra_list$Pongo_pygmaeus
      )
    return(spectra)
  }

normalize_spectra =
  function(spectra, trinuc_composition){
    spectra_norm =
      sapply(
        colnames(
          spectra
        ),
        function(x)
          (spectra[, x] /
             rowSums(spectra)) /
          (trinuc_composition[[substr(x, 1, 3)]] /
             sum(trinuc_composition))
      )
    return(spectra_norm)
  }

reweight_indiv_spectra <- 
  function(spectra, trinuc_composition, reweight_trinuc_composition, fraction=TRUE){
    species_df_reweight = 
      sweep(spectra, 
            2, 
            unlist((reweight_trinuc_composition/trinuc_composition)[substr(colnames(spectra), 1, 3)]),
            `*`)
    if(fraction){
      return(species_df_reweight / rowSums(species_df_reweight))
    }else{
      return(species_df_reweight)
    }
  }

reweight_indiv_spectra_for_chisq <-
  function(spectra, trinuc_composition, alt_trinuc_composition){
    spectra_reweight = 
      sapply(
        collapsed_trinuc_mutations,
        function(m)
          if(trinuc_composition[substr(m, 1, 3)] > 
             alt_trinuc_composition[substr(m, 1, 3)]){
            spectra[, m] * 
              (alt_trinuc_composition[[substr(m, 1, 3)]] /
                 trinuc_composition[[substr(m, 1, 3)]])
          } else {
            spectra[, m]
          }
      )
    rownames(spectra_reweight) = rownames(spectra)
    return(spectra_reweight)
  }

combi <- function(vec1)
{
  si <- length(vec1)
  first <- rep(vec1, (si-1):0)
  secR <- rev(vec1)
  second <- secR[sequence(1:(si-1))]
  second <- rev(second)
  combi <- matrix(cbind(first, second), ncol = 2)
  return(combi)
}

compare_two_species_indiv_spectra_chisq <-
  function(c1_spectra, c1_species, c2_spectra, c2_species){
    if(all(c1_spectra == c2_spectra) & c1_species == c2_species){
        indiv_list = which(startsWith(rownames(c1_spectra), c1_species))
        indiv_combi = combi(indiv_list)
      } else {
        indiv_list_c1 = which(startsWith(rownames(c1_spectra), c1_species))
        indiv_list_c2 = which(startsWith(rownames(c2_spectra), c2_species))
        indiv_combi = expand.grid(indiv_list_c1, indiv_list_c2)
      }
    return(
          apply(
            indiv_combi,
            1,
            function(x)
              # print(unlist(c1_spectra[x[1],]))
              chisq.test(
                matrix(
                  c(
                    unlist(c1_spectra[x[1],]),
                    unlist(c2_spectra[x[2],])
                  ),
                  nrow=2, byrow = TRUE
                )
              )$statistic
          )
      )
  }

species_initials <- c("H", "PanT", "PanP", "G", "PonA", "PonP")
names(species_initials) <- species

combi_species <- 
  rbind(
    matrix(c(species, species), ncol = 2),
    combi(species)
  )
combi_species_levels <-
  c("G|PonP", "PanP|PonP", "PanT|PonP", "H|PonP", 
    "G|PonA", "PanP|PonA", "PanT|PonA", "H|PonA", 
    "PanP|G", "PanT|G", "H|G", 
    "H|PanP", "H|PanT",
    "PanT|PanP",
    "PonA|PonP",
    "PonP|PonP", "PonA|PonA", "G|G", "PanP|PanP", "PanT|PanT", "H|H")

compare_two_indiv_spectra_chisq <-
  function(c1_spectra, c1_trinuc_composition, c2_spectra, c2_trinuc_composition){
    c1_spectra_reweight = 
      reweight_indiv_spectra_for_chisq(
        c1_spectra, c1_trinuc_composition, c2_trinuc_composition
      )
    c2_spectra_reweight = 
      reweight_indiv_spectra_for_chisq(
        c2_spectra, c2_trinuc_composition, c1_trinuc_composition
      )
    chi2_list_of_df = 
      apply(
        combi_species,
        1,
        function(x){
          chi2_stats = 
            compare_two_species_indiv_spectra_chisq(
              c1_spectra_reweight, x[1], c2_spectra_reweight, x[2]
            )
          return(
            data.frame(
              s1 = rep(x[1], length(chi2_stats)),
              s2 = rep(x[2], length(chi2_stats)),
              chi2 = chi2_stats
            )
          )
        }
      )
    chisq_df = ldply(chi2_list_of_df, data.frame)
    chisq_df$"comparison" = 
      factor(
        paste(
          species_initials[chisq_df$s1], 
          species_initials[chisq_df$s2], 
          sep = "|"),
        levels = combi_species_levels
        )
    return(
      chisq_df
    )
  }

plot_indiv_compartment_chisq <-
  function(c1_spectra, c1_trinuc_composition, c2_spectra, c2_trinuc_composition, figure_title, species_list){
    chisq_df = compare_two_indiv_spectra_chisq(
      c1_spectra, c1_trinuc_composition,
      c2_spectra, c2_trinuc_composition
    )
    ggplot(
      subset(chisq_df, chisq_df$s1 %in% species_list & chisq_df$s2 %in% species_list), 
      aes(x=log(chi2), y=comparison, group=comparison)) +
      geom_density_ridges(scale = 10, size = 0.25, rel_min_height = 0.03) +
      theme_ridges() +
      labs(title = figure_title) +
      theme(axis.text.x = element_text(angle=90))
  }

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Using opacity to show subspecies - this is tough
# species_colors <- gg_color_hue(6)
# subspecies_colors_rgb <- species_colors[c(1, 1, 1, 2, 3, 4, 4, 4, 4, 5, 6)]
# subspecies_colors <- 
#   paste0(subspecies_colors_rgb, 
#          as.hexmode(c(85, 170, 255, 255, 255, 63, 127, 191, 255, 255, 255)))

# Using saturation and brightness - this is tougher
species_colors <- c("#F77167", "#9B8600", "#00BA38", "#0095A3", "#619CFF", "#F564E3")
subspecies_colors <- c("#F99A93", "#F77167", "#A54F48", "#9B8600", "#00BA38", "#00E5ED", "#00C5CC", "#0095A3", "#00597A", "#619CFF", "#F564E3")



load_species_3mer_spectrum <-
  function(compartment_name, nmer_dir = "./snp_data/nmer_mutation_counts/"){
    species_spectra_df = 
      data.frame(
        sapply(species, function(s){
          filename = paste0(nmer_dir, s, "_", compartment_name, "_species_3mer_counts.txt")
          return(unlist(read.table(filename, header = T)))
        })
      )
    row.names(species_spectra_df) = collapsed_trinuc_mutations
    return(species_spectra_df)
  }

normalize_species_spectra <-
  function(spectra, trinuc_composition){
    species_df_norm = 
      as.data.frame(
        sapply(
          colnames(spectra),
          function(s)
            (spectra[,c(s)] / 
               sum(spectra[,s])) /
            (trinuc_composition[substr(collapsed_trinuc_mutations, 1, 3)] /
               sum(trinuc_composition))
        )
      )
    rownames(species_df_norm) <- collapsed_trinuc_mutations
    return(species_df_norm)
  }

reweight_species_spectra <-
  function(spectra, trinuc_composition, reweight_trinuc_composition, fraction=TRUE, for_chi_sq=FALSE){
    species_df_reweight = 
      as.data.frame(
        sapply(
          species,
          function(s){
            if(for_chi_sq){
              return(
                unlist(
                  sapply(
                    collapsed_trinuc_mutations,
                    function(m)
                      if(
                        reweight_trinuc_composition[substr(m, 1, 3)] / 
                        trinuc_composition[substr(m, 1, 3)] < 1
                      ){
                        return(
                          as.numeric(spectra[m, s]) * 
                            reweight_trinuc_composition[substr(m, 1, 3)] / 
                            trinuc_composition[substr(m, 1, 3)])
                      }else{
                        return(spectra[m, s])
                      }
                  )
                )
              )
            }
            else{
              reweighted_counts = 
                unlist(
                  spectra[,c(s)] *
                    (reweight_trinuc_composition[substr(collapsed_trinuc_mutations, 1, 3)] /
                       trinuc_composition[substr(collapsed_trinuc_mutations, 1, 3)])
                )
              if(fraction){
                return(reweighted_counts / sum(reweighted_counts))
              }
              else{
                return(reweighted_counts)
              }
            }
          }
        )
      )
    rownames(species_df_reweight) <- collapsed_trinuc_mutations
    return(species_df_reweight)
  }


chi_square_species_spectra <-
  function(c1_spectra, c1_trinuc_composition, c2_spectra, c2_trinuc_composition){
    normalized_counts1 = 
      reweight_species_spectra(
        c1_spectra,
        c1_trinuc_composition,
        c2_trinuc_composition,
        for_chi_sq = T
      )
    normalized_counts2 = 
      reweight_species_spectra(
        c2_spectra,
        c2_trinuc_composition,
        c1_trinuc_composition,
        for_chi_sq = T
      )
    chi_square_df = 
      sapply(
        species,
        function(s)
          sapply(
            collapsed_trinuc_mutations,
            function(m)
              chisq.test(
                matrix(
                  c(
                    normalized_counts1[m, s],
                    normalized_counts2[m, s],
                    sum(normalized_counts1[, s]),
                    sum(normalized_counts2[, s])
                  ), nrow = 2, ncol = 2
                )
              )$p.value
          )
      )
    return(chi_square_df)
  }

log_ratio_spectra <-
  function(spectra1, spectra2){
    log_ratio_df = 
      as.data.frame(
        sapply(
          species,
          function(s)
            log(
              unlist(spectra1[, s]) /
                unlist(spectra2[, s])
            )
        )
      )
    rownames(log_ratio_df) = collapsed_trinuc_mutations
    return(log_ratio_df)
  }

log_ratio_spectra <-
  function(spectra1, spectra2){
    log_ratio_df = 
      as.data.frame(
        sapply(
          species,
          function(s)
            log(
              unlist(spectra1[, s]) /
                unlist(spectra2[, s])
            )
        )
      )
    rownames(log_ratio_df) = collapsed_trinuc_mutations
    return(log_ratio_df)
  }

log_ratio_species <-
  function(spectra){
    log_ratio_df =
      as.data.frame(
        apply(
          combi(species),
          1,
          function(s)
            log(
              unlist(spectra[, s[1]]) /
                unlist(spectra[, s[2]])
            )
        )
      )
    rownames(log_ratio_df) = collapsed_trinuc_mutations
    colnames(log_ratio_df) = paste(species_initials[combi(species)[,1]],
                                   species_initials[combi(species)[,2]],
                                   sep = "|")
    return(log_ratio_df)
  }

generate_heatmap_plot_single_species <- function(spectra1, spectra2, s) {
  db = log_ratio_spectra(spectra1, spectra2)
  db$"five_prime_and_mut" = five_prime_and_mut
  db$"three_prime" = three_prime
  plot_list =
    lapply(
      rev(unique(substr(levels(db[, c("five_prime_and_mut")]), 2, 4))),
      function(mut)
        ggplot(subset(db, substr(db[, c("five_prime_and_mut")], 2, 4) == mut), aes(three_prime, five_prime_and_mut)) +
        geom_tile(aes(fill=subset(db[, c(s)], substr(db[, c("five_prime_and_mut")], 2, 4) == mut)), color = "white") +
        scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = 0, limits = c(min(db[,s]), max(db[,s]))) +
        theme_minimal() +
        coord_equal() +
        ylab(sub("\\.", ">", mut)) +
        scale_y_discrete(labels=c("A", "C", "G", "T")) +
        theme(axis.text.x = element_blank(), 
              axis.ticks = element_blank(), 
              axis.title.x = element_blank(), 
              axis.title.y = element_text(angle = 0, vjust = 0.5),
              legend.position = "none",
              panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "white"),
              plot.margin =unit(c(0,0,0,0), "cm"))
    )
  plot_list[[1]] =
    plot_list[[1]] +
    ggtitle(s) +
    theme(plot.title=element_text())
  plot_list[[6]] = 
    plot_list[[6]] + 
    theme(axis.text.x=element_text())
  for_legend = 
    ggplot(db, aes(three_prime, five_prime_and_mut)) +
    geom_tile(aes(fill=db[, c(s)]), color = "white") +
    scale_fill_gradient2(low = muted("blue"), 
                         mid = "white", 
                         high = muted("red"), 
                         midpoint = 0, 
                         limits = c(min(db[,s]), max(db[,s]))) +
    theme(legend.title = element_blank())
  plot_list[[7]] <- cowplot::get_legend(for_legend)
  hlay = rbind(c(1,7),
               c(2,7),
               c(3,7),
               c(4,7),
               c(5,7),
               c(6,7))
  grid.arrange(grobs=plot_list, 
               ncol = 2, 
               layout_matrix = hlay,
               heights = unit(c(2.6, rep(2, 4), 2.25), "cm"))
}

generate_heatmap_plot_multiple_species <- function(spectra1, spectra2, species_list=species) {
  db = log_ratio_spectra(spectra1, spectra2)
  db$"five_prime_and_mut" = five_prime_and_mut
  db$"three_prime" = three_prime
  plot_list_of_lists =
    lapply(
      species_list,
      function(s){
        species_plot_list = 
          lapply(
            rev(unique(substr(levels(db[, c("five_prime_and_mut")]), 2, 4))),
            function(mut)
              ggplot(subset(db, substr(db[, c("five_prime_and_mut")], 2, 4) == mut), aes(three_prime, five_prime_and_mut)) +
              geom_tile(aes(fill=subset(db[, c(s)], substr(db[, c("five_prime_and_mut")], 2, 4) == mut)), color = "white") +
              scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = 0, limits = c(min(db[,species_list]), max(db[,species_list]))) +
              theme_minimal() +
              coord_equal() +
              # ylab(sub("\\.", ">", mut)) +
              # scale_y_discrete(labels=c("A", "C", "G", "T")) +
              theme(axis.text.x = element_blank(), 
                    axis.text.y = element_blank(),
                    axis.ticks = element_blank(), 
                    axis.title.x = element_blank(), 
                    axis.title.y = element_blank(),
                    # axis.title.y = element_text(angle = 0, vjust = 0.5),
                    legend.position = "none",
                    panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "white"),
                    plot.margin =unit(c(0,0,0,0), "cm"))
          )
        species_plot_list[[1]] =
          species_plot_list[[1]] +
          ggtitle(s) +
          theme(plot.title=element_text())
        species_plot_list[[6]] = 
          species_plot_list[[6]] + 
          theme(axis.text.x=element_text())
        return(species_plot_list)
      }
    )
  for_legend = 
    ggplot(db, aes(three_prime, five_prime_and_mut)) +
    geom_tile(aes(fill=db[, c(species_list[1])]), color = "white") +
    scale_fill_gradient2(low = muted("blue"), 
                         mid = "white", 
                         high = muted("red"), 
                         midpoint = 0, 
                         limits = c(min(db[,species_list]), max(db[,species_list]))) +
    theme(legend.title = element_blank())
  plot_list = unlist(plot_list_of_lists, recursive = F)
  plot_list[[6*length(species_list) + 1]] = cowplot::get_legend(for_legend)
  mut_levels = rev(unique(substr(levels(db[, c("five_prime_and_mut")]), 2, 4)))
  # Adding back the y labels for the farthest left plot
  for(i in c(1:length(mut_levels))){
    plot_list[[i]] = 
      plot_list[[i]] + 
      ylab(sub("\\.", ">", mut_levels[i])) +
      scale_y_discrete(labels=c("A", "C", "G", "T")) +
      theme(axis.text.y=element_text(),
            axis.title.y = element_text(angle = 0, vjust = 0.5))
  }
  hlay = cbind(matrix(c(1:(6*length(species_list))), nrow=6), 6*length(species_list) + 1)
  grid.arrange(grobs=plot_list, 
               ncol = 7, 
               layout_matrix = hlay,
               heights = unit(c(2.6, rep(2, 4), 2.25), "cm"))
}

generate_heatmap_plot_single_species_chisq <- function(c1_spectra, c1_trinuc_composition, c2_spectra, c2_trinuc_composition, s) {
  c1_spectra_reweight = reweight_species_spectra(c1_spectra, c1_trinuc_composition, c2_trinuc_composition, for_chi_sq = TRUE)
  c2_spectra_reweight = reweight_species_spectra(c2_spectra, c2_trinuc_composition, c1_trinuc_composition, for_chi_sq = TRUE)
  c1_spectra_reweight_rates = sweep(c1_spectra_reweight, 2, colSums(c1_spectra_reweight), `/`)
  c2_spectra_reweight_rates = sweep(c2_spectra_reweight, 2, colSums(c2_spectra_reweight), `/`)
  db = log_ratio_spectra(c1_spectra_reweight_rates, c2_spectra_reweight_rates)
  db$"five_prime_and_mut" = five_prime_and_mut
  db$"three_prime" = three_prime
  chisq_db = chi_square_species_spectra(c1_spectra, c1_trinuc_composition, c2_spectra, c2_trinuc_composition) < 1e-4
  colnames(chisq_db) = paste0(colnames(chisq_db), "_chisq")
  db = cbind(db, chisq_db)
  plot_list =
    lapply(
      rev(unique(substr(levels(db[, c("five_prime_and_mut")]), 2, 4))),
      function(mut)
        ggplot(subset(db, substr(db[, c("five_prime_and_mut")], 2, 4) == mut), aes(three_prime, five_prime_and_mut)) +
        geom_tile(aes(fill=subset(db[, c(s)], substr(db[, c("five_prime_and_mut")], 2, 4) == mut)), color = "white") +
        scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = 0, limits = c(min(db[,s]), max(db[,s]))) +
        theme_minimal() +
        coord_equal() +
        ylab(sub("\\.", ">", mut)) +
        scale_y_discrete(labels=c("A", "C", "G", "T")) +
        theme(axis.text.x = element_blank(), 
              axis.ticks = element_blank(), 
              axis.title.x = element_blank(), 
              axis.title.y = element_text(angle = 0, vjust = 0.5),
              legend.position = "none",
              panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "white"),
              plot.margin =unit(c(0,0,0,0), "cm")) +
        geom_point(data = subset(db, substr(db[, c("five_prime_and_mut")], 2, 4) == mut & 
                                     db[, paste0(s, "_chisq")]),
                   aes(x=three_prime, y=five_prime_and_mut), size = 2)
    )
  plot_list[[1]] =
    plot_list[[1]] +
    ggtitle(s) +
    theme(plot.title=element_text())
  plot_list[[6]] = 
    plot_list[[6]] +  
    theme(axis.text.x=element_text())
  for_legend = 
    ggplot(db, aes(three_prime, five_prime_and_mut)) +
    geom_tile(aes(fill=db[, c(s)]), color = "white") +
    scale_fill_gradient2(low = muted("blue"), 
                         mid = "white", 
                         high = muted("red"), 
                         midpoint = 0, 
                         limits = c(min(db[,s]), max(db[,s]))) +
    theme(legend.title = element_blank())
  plot_list[[7]] <- cowplot::get_legend(for_legend)
  hlay = rbind(c(1,7),
               c(2,7),
               c(3,7),
               c(4,7),
               c(5,7),
               c(6,7))
  grid.arrange(grobs=plot_list, 
               ncol = 2, 
               layout_matrix = hlay,
               heights = unit(c(2.6, rep(2, 4), 2.25), "cm"))
}

# Spectra need to be reweighted
log_odds_spectra_correlation_species_heatmap <- function(c1_spectra, c2_spectra){
  lr_df = log_ratio_spectra(c1_spectra, c2_spectra)
  pval_matrix = 
    sapply(
      species,
      function(s1)
        sapply(
          species,
          function(s2)
            cor.test(lr_df[,s1], lr_df[,s2])$p.value
        )
    )
  pval_df = 
    data.frame(
      s1 = factor(rep(species, each = 6), levels=c("Pongo_abelii", "Pongo_pygmaeus", "Gorilla", "Homo", "Pan_troglodytes", "Pan_paniscus")),
      s2 = factor(rep(species, times = 6), levels=c("Pan_paniscus", "Pan_troglodytes", "Homo", "Gorilla", "Pongo_pygmaeus", "Pongo_abelii")),
      pval = as.vector(pval_matrix)
    )
  return(
    ggplot(pval_df, aes(s1, s2)) +
      geom_tile(aes(fill=-log(pval)), color = "white") +
      scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = 0) +
      theme_minimal() +
      scale_x_discrete(position = "top") +
      coord_equal()
  )
}

# Spectra need to be reweighted
log_odds_spectra_correlation_species <- function(c1_spectra, c2_spectra){
  db = log_ratio_spectra(c1_spectra, c2_spectra)
  corr_db = 
    data.frame(
      comparison = paste(species_initials[combi(species)[,1]],
                         species_initials[combi(species)[,2]],
                         sep = "|"),
      rho = apply(
        combi(species),
        1,
        function(s){
          cor.test(db[,s[1]], db[,s[2]])$estimate
        }
      ),
      pval = apply(
        combi(species),
        1,
        function(s){
          cor.test(db[,s[1]], db[,s[2]])$p.value
        }
      )
    )
  return(
    corr_db
  )
}

log_odds_species_correlation_spectra_heatmap <- function(c1_spectra, c2_spectra){
  c1_lr = log_ratio_species(c1_spectra)
  c2_lr = log_ratio_species(c2_spectra)
  pval_list = 
    sapply(
      paste(species_initials[combi(species)[,1]],
            species_initials[combi(species)[,2]],
            sep = "|"),
      function(s){
        cor.test(c1_lr[,s], c2_lr[,s])$p.value
      }
    )
  pval_df = 
    data.frame(
      s1=factor(sub("\\|.*", "", names(pval_list)), levels=c("PonP", "PonA", "G", "H", "PanT", "PanP")),
      s2=factor(sub("^.*\\|", "", names(pval_list)), levels=c("PanP", "PanT", "H", "G", "PonA", "PonP")),
      pval = pval_list
    )
  return(
    ggplot(pval_df, aes(s1, s2)) +
      geom_tile(aes(fill=-log(pval)), color = "white") +
      scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = 0) +
      theme_minimal() +
      scale_x_discrete(position = "top") +
      coord_equal()
  )
}

log_odds_species_correlation_spectra <- function(c1_spectra, c2_spectra){
  c1_lr = log_ratio_species(c1_spectra)
  c2_lr = log_ratio_species(c2_spectra)
  return(
    sapply(
      paste(species_initials[combi(species)[,1]],
            species_initials[combi(species)[,2]],
            sep = "|"),
      function(s){
        cor.test(c1_lr[,s], c2_lr[,s])$p.value
      }
    )
  )
}

log_odds_species_correlation_spectra_table <- function(c1_spectra, c2_spectra){
  c1_lr = log_ratio_species(c1_spectra)
  c2_lr = log_ratio_species(c2_spectra)
  return(
    sapply(
      paste(species_initials[combi(species)[,1]],
            species_initials[combi(species)[,2]],
            sep = "|"),
      function(s){
        cor.test(c1_lr[,s], c2_lr[,s])$p.value
      }
    )
  )
}

generate_pc_axis_label <- function(pca, pc_axis){
  axis_label = 
    paste0(
      "PC", 
      pc_axis, 
      " (",
      round(summary(pca)$importance[2,pc_axis] * 100, digits = 1),
      "%)"
    )
  return(axis_label)
}


