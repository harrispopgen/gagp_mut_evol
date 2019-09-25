# Code for analysis and extraction of ERVs from great ape sequences
# 2018-07-16

options(stringsAsFactors = F)
require(data.table)

setwd("/net/harris/vol1/project/primate_ervs")
# setwd("~/Documents/Harris_project_directory")

rm_header <- read.table("./repeatmasker_header.txt")[,1]
# gorilla_data <- fread("./erv_data/gorGor5_ReapeatMasker.txt", verbose=TRUE)
# human_data <- fread("./erv_data/hg38_ReapeatMasker.txt", verbose=TRUE)
# bonobo_data <- fread("./erv_data/panPan2_repeatMasker.txt", verbose=TRUE)
hg18_data <- fread("/net/harris/vol1/data/hg18/hg18_repeatMasker.txt", verbose=TRUE)

# names(gorilla_data) <- rm_header
# names(human_data) <- rm_header
# names(bonobo_data) <- rm_header
names(hg18_data) <- rm_header

# gorilla_data$N <- row.names(gorilla_data)
# human_data$N <- row.names(human_data)
# bonobo_data$N <- row.names(bonobo_data)
hg18_data$N <- row.names(hg18_data)

# gorilla_ltr_data <- subset(gorilla_data, gorilla_data$repClass == "LTR")
# human_ltr_data <- subset(human_data, human_data$repClass == "LTR")
# bonobo_ltr_data <- subset(bonobo_data, bonobo_data$repClass == "LTR")
hg18_ltr_data <- subset(hg18_data, hg18_data$repClass == "LTR")

erv_family_subset <- c("ERV1", "ERVK", "ERVL", "Gypsy")

# gorilla_all_erv <- subset(gorilla_ltr_data, gorilla_ltr_data$repFamily %in% erv_family_subset)
# gorilla_malr <- subset(gorilla_ltr_data, gorilla_ltr_data$repFamily == "ERVL-MaLR")
# gorilla_all_erv_bed <- gorilla_all_erv[, c("genoName", "genoStart", "genoEnd", "N")]
# gorilla_malr_bed <- gorilla_malr[, c("genoName", "genoStart", "genoEnd", "N")]

# human_all_erv <- subset(human_ltr_data, human_ltr_data$repFamily %in% erv_family_subset)
# human_malr <- subset(human_ltr_data, human_ltr_data$repFamily %in% "ERVL-MaLR")
# human_all_erv_bed <- human_all_erv[, c("genoName", "genoStart", "genoEnd", "N")]
# human_malr_bed <- human_malr[, c("genoName", "genoStart", "genoEnd", "N")]

# bonobo_all_erv <- subset(bonobo_ltr_data, bonobo_ltr_data$repFamily %in% erv_family_subset)
# bonobo_malr <- subset(bonobo_ltr_data, bonobo_ltr_data$repFamily %in% "ERVL-MaLR")
# bonobo_all_erv_bed <- bonobo_all_erv[, c("genoName", "genoStart", "genoEnd", "N")]
# bonobo_malr_bed <- bonobo_malr[, c("genoName", "genoStart", "genoEnd", "N")]

hg18_all_erv <- subset(hg18_ltr_data, hg18_ltr_data$repFamily %in% erv_family_subset)
hg18_malr <- subset(hg18_ltr_data, hg18_ltr_data$repFamily %in% "MaLR") # Nomenclature is different for hg18 I think?
hg18_all_erv_bed <- hg18_all_erv[, c("genoName", "genoStart", "genoEnd", "N")]
hg18_malr_bed <- hg18_malr[, c("genoName", "genoStart", "genoEnd", "N")]


# write.table(gorilla_all_erv, "./erv_data/gorGor_all_erv.txt", sep = "\t")
# write.table(gorilla_malr, "./erv_data/gorGor_malr.txt", sep = "\t")
# write.table(gorilla_all_erv_bed, "./erv_data/gorGor_all_erv.bed", sep = "\t", row.names = F, quote = F, col.names = F)
# write.table(gorilla_malr_bed, "./erv_data/gorGor_malr.bed", sep = "\t", row.names = F, quote = F, col.names = F)
# write.table(human_all_erv, "./erv_data/hg38_all_erv.txt", sep = "\t")
# write.table(human_malr, "./erv_data/hg38_malr.txt", sep = "\t")
# write.table(human_all_erv_bed, "./erv_data/hg38_all_erv.bed", sep = "\t", row.names = F, quote = F, col.names = F)
# write.table(human_malr_bed, "./erv_data/hg38_malr.bed", sep = "\t", row.names = F, quote = F, col.names = F)
# write.table(bonobo_all_erv, "./erv_data/panPan2_all_erv.txt", sep = "\t")
# write.table(bonobo_malr, "./erv_data/panPan2_malr.txt", sep = "\t")
# write.table(bonobo_all_erv_bed, "./erv_data/panPan2_all_erv.bed", sep = "\t", row.names = F, quote = F, col.names = F)
# write.table(bonobo_malr_bed, "./erv_data/panPan2_malr.bed", sep = "\t", row.names = F, quote = F, col.names = F)
write.table(hg18_all_erv, "./erv_data/hg18_all_erv.txt", sep = "\t")
write.table(hg18_malr, "./erv_data/hg18_malr.txt", sep = "\t")
write.table(hg18_all_erv_bed, "./erv_data/hg18_all_erv.bed", sep = "\t", row.names = F, quote = F, col.names = F)
write.table(hg18_malr_bed, "./erv_data/hg18_malr.bed", sep = "\t", row.names = F, quote = F, col.names = F)


