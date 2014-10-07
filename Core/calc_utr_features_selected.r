#! /usr/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)

headersTXT = c(
"Miranda_score",
"miR_ID",
"mRNA_ID",
"Start_position",
"End_position",
"Seed_match_6mer2",
"miR_match_P01",
"Seed_match_7mer2",
"Seed_match_7mer1",
"Seed_MFE",
"X3p_MFE",
"Target_UC_comp",
"miR_match_P09",
"miR_match_P02",
"Seed_GU",
"miR_match_P07",
"miR_match_P19",
"miR_match_P15"
)

read.table(args[1],sep="\t", header=TRUE) -> features
features = subset(features, select = headersTXT)

#library(seqinr)
#read.fasta(file=args[2], as.string=TRUE) -> fa
#fa=fa[which(!duplicated(names(fa)))]
#utr_length = t(sapply(names(fa), function(i) {
#	c(i, nchar(fa[[i]][1]))
#}))

#remove records with UTRs not in fasta file
#features = subset(features, mRNA_ID %in% names(fa))

attach(features)
features_numeric = subset(features, select = -c(miR_ID,mRNA_ID))
aggregate(features_numeric, by=list(miR_ID,mRNA_ID), sum) -> features_sum
aggregate(1:dim(features_numeric)[1], by=list(miR_ID,mRNA_ID), length) -> number_sites
colnames(number_sites)[3] <- "number_sites"
aggregate(features_numeric, by=list(miR_ID,mRNA_ID), mean) -> features_mean
aggregate(features_numeric, by=list(miR_ID,mRNA_ID), max) -> features_max
aggregate(features_numeric, by=list(miR_ID,mRNA_ID), min) -> features_min
detach(features)

features_sum2 = features_sum
features_mean2 = features_mean
features_max2 = features_max
features_min2 = features_min

colnames(features_sum2) <- paste0(colnames(features_sum2), ".sum")
colnames(features_mean2) <- paste0(colnames(features_mean2), ".mean")
colnames(features_max2) <- paste0(colnames(features_max2), ".max")
colnames(features_min2) <- paste0(colnames(features_min2), ".min")
output = merge(number_sites, features_sum2, by=c(1,2))
output = merge(output, features_mean2, by=c(1,2))
output = merge(output, features_max2, by=c(1,2))
output = merge(output, features_min2, by=c(1,2))

selected_features = 
c("Group.1",
"Group.2",
"Miranda_score.max",
"Seed_match_6mer2.mean",
"miR_match_P01.min",
"Seed_match_7mer2.max",
"Seed_match_7mer1.mean",
"Seed_MFE.min",
"X3p_MFE.mean",
"Target_UC_comp.mean",
"miR_match_P09.mean",
"miR_match_P02.min",
"Seed_GU.mean",
"miR_match_P07.mean",
"Start_position.min",
"miR_match_P19.min",
"miR_match_P15.min"
)

output2 <- output[,selected_features]
write.table(output2, file=paste0(args[2], ".csv"), sep=",", row.names=FALSE)
