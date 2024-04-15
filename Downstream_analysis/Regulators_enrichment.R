library(data.table)
library(dplyr)

#RBP
for(dir in c("m6A2me", "me2m6A")){
	RBP_enrichment <- fread(paste0("~/Downstream_analysis/Regulators/Lung.m6A_peak_RBP.", dir, ".txt"), sep = "\t", header = FALSE, col.names = c("RBP", "signif", "background")) %>%
		mutate(prop_signif = signif/sum(signif), prop_background = background/sum(background), Fold = prop_signif/prop_background) %>%
		mutate(hyper.p.right = phyper(signif - 1, background, sum(background) - background, sum(signif), lower.tail = FALSE)) %>%
		mutate(fdr.hyper.right = p.adjust(hyper.p.right, method = "fdr", n = length(hyper.p.right)))
	#fisher test
	RBP_enrichment$fisher.p.greater <- 0
	for(i in 1:nrow(RBP_enrichment)){
		signif <- as.numeric(RBP_enrichment[i, "signif"])
		background <- as.numeric(RBP_enrichment[i, "background"])
		RBP_enrichment[i, "fisher.p.greater"] <- fisher.test(data.frame(c1 = c(signif, background - signif),
										c2 = c(sum(RBP_enrichment$signif) - signif, sum(RBP_enrichment$background) - background - (sum(RBP_enrichment$signif) - signif))), alternative = "greater")[["p.value"]]
	}
	RBP_enrichment <- RBP_enrichment %>% mutate(fdr.fisher.greater = p.adjust(fisher.p.greater, method = "fdr", n = length(fisher.p.greater)))

	write.table(RBP_enrichment, file = paste0("~/Downstream_analysis/Regulators/Lung.m6A_peak_RBP.", dir, ".pval.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
	write.table(RBP_enrichment %>% filter(fdr.fisher.greater < 0.05), file = paste0("~/Downstream_analysis/Regulators/Lung.m6A_peak_RBP.", dir, ".signif.list"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

#TF, integrate chromosome split results
for(dir in c("m6A2me", "me2m6A")){
	TF_enrichment <- data.frame()
	for(i in 1:22){
		TF_enrichment0 <- fread(paste0("~/Downstream_analysis/Regulators/Lung_DNAme_split_bed/Lung.DNAme_TF.", dir, "_chr", i, ".txt"), sep = "\t", header = FALSE, col.names = c("TF", "signif_chr", "background_chr"))
		TF_enrichment <- bind_rows(TF_enrichment, TF_enrichment0)
	}
	TF_enrichment <- TF_enrichment %>% group_by(TF) %>% summarise(signif = sum(signif_chr), background = sum(background_chr)) %>%
		mutate(prop_signif = signif/sum(signif), prop_background = background/sum(background), Fold = prop_signif/prop_background) %>%
		mutate(hyper.p.right = phyper(signif - 1, background, sum(background) - background, sum(signif), lower.tail = FALSE)) %>%
		mutate(fdr.hyper.right = p.adjust(hyper.p.right, method = "fdr", n = length(hyper.p.right)))
	##fisher test
	TF_enrichment$fisher.p.greater <- 0
	for(i in 1:nrow(TF_enrichment)){
		signif <- as.numeric(TF_enrichment[i, "signif"])
                background <- as.numeric(TF_enrichment[i, "background"])
		TF_enrichment[i, "fisher.p.greater"] <- fisher.test(data.frame(c1 = c(signif, background - signif),
									       c2 = c(sum(TF_enrichment$signif) - signif, sum(TF_enrichment$background) - background - (sum(TF_enrichment$signif) - signif))), alternative = "greater")[["p.value"]]
	}
	TF_enrichment <- TF_enrichment %>% mutate(fdr.fisher.greater = p.adjust(fisher.p.greater, method = "fdr", n = length(fisher.p.greater)))

	write.table(TF_enrichment, file = paste0("~/Downstream_analysis/Regulators/Lung.DNAme_TF.", dir, ".pval.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
	write.table(TF_enrichment %>% filter(fdr.fisher.greater < 0.05), file = paste0("~/Downstream_analysis/Regulators/Lung.DNAme_TF.", dir, ".signif.list"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}
