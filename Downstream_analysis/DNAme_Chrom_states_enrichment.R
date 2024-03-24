library(data.table)
library(dplyr)

for(tis in c("Lung", "Muscle")){
	for(dir in c("m6A2me", "me2m6A")){
		DNAme_enrichment <- fread(paste0("~/Downstream_analysis/Chrom_states_enrichment/", tis, ".DNAme.", dir, ".txt"), sep = "\t", header = FALSE, col.names = c("Chrom_states", "signif", "background")) %>%
			mutate(prop_signif = signif/sum(signif), prop_background = background/sum(background), Fold = prop_signif/prop_background) %>%
			mutate(hyper.p.left = phyper(signif, background, sum(background) - background, sum(signif), lower.tail = TRUE),
			       hyper.p.right = phyper(signif - 1, background, sum(background) - background, sum(signif), lower.tail = FALSE))
		#fisher test
		DNAme_enrichment$$fisher.p.twosided <- 0
		for(i in 1:nrow(DNAme_enrichment)){
			signif <- as.numeric(DNAme_enrichment[i, "signif"])
			background <- as.numeric(DNAme_enrichment[i, "background"])
			DNAme_enrichment[i, "fisher.p.twosided"] <- fisher.test(data.frame(c1 = c(signif, background - signif),
											   c2 = c(sum(DNAme_enrichment$signif) - signif, sum(DNAme_enrichment$background) - background - (sum(DNAme_enrichment$signif) - signif))),
							      			alternative = "two.sided")[["p.value"]]
		}
		DNAme_enrichment <- DNAme_enrichment %>% mutate(fdr.fisher.twosided = p.adjust(fisher.p.twosided, method = "fdr", n = length(fisher.p.twosided)))

		write.table(DNAme_enrichment, file = paste0("~/Downstream_analysis/Chrom_states_enrichment/", tis, ".DNAme.", dir, ".pval.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
	}
}
