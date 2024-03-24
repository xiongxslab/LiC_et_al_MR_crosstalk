library(data.table)
library(dplyr)

for(tis in c("Lung", "Muscle")){
	m6A_anno <- fread(paste0("~/Downstream_analysis/m6A_Genomic_enrichment/", tis, ".m6A_peaks.merged_tested.anno.txt"), sep = "\t", header = TRUE) %>%
		rename(Peak_ID = matches("PeakID")) %>%
		select(Peak_ID, Annotation) %>%
		mutate(Annotation = gsub(" \\(.*", "", Annotation)) %>%
		filter(Annotation != "Intergenic" & Annotation != "intron")

	for(dir in c("m6A2me", "me2m6A")){
		m6A_signif <- fread(paste0("~/2SMR/Results/Lung.m6A_peak.", dir, ".signif.list"), sep = "\t", header = FALSE, col.names = c("Peak_ID"))
		m6A_background <- fread(paste0("~/2SMR/Results/Lung.m6A_peak.", dir, ".tested.list"), sep = "\t", header = FALSE, col.names = c("Peak_ID"))

		m6A_enrichment <- right_join(m6A_anno %>% inner_join(y = m6A_signif, by = "Peak_ID") %>% group_by(Annotation) %>% summarise(signif = n()),
					     m6A_anno %>% inner_join(y = m6A_background, by = "Peak_ID") %>% group_by(Annotation) %>% summarise(background = n()), by = "Annotation") %>% replace_na(list(signif = 0))

		m6A_enrichment <- m6A_enrichment %>%
			mutate(prop_signif = signif/sum(signif), prop_background = background/sum(background), Fold = prop_signif/prop_background) %>%
			mutate(hyper.p.left = phyper(signif, background, sum(background) - background, sum(signif), lower.tail = TRUE),
			       hyper.p.right = phyper(signif - 1, background, sum(background) - background, sum(signif), lower.tail = FALSE))

		#fisher test
		m6A_enrichment$$fisher.p.twosided <- 0
		for(i in 1:nrow(m6A_enrichment)){
			signif <- as.numeric(m6A_enrichment[i, "signif"])
                        background <- as.numeric(m6A_enrichment[i, "background"])
                        m6A_enrichment[i, "fisher.p.twosided"] <- fisher.test(data.frame(c1 = c(signif, background - signif),
											 c2 = c(sum(m6A_enrichment$signif) - signif, sum(m6A_enrichment$background) - background - (sum(m6A_enrichment$signif) - signif))), alternative = "two.sided")[["p.value"]]
                }
                m6A_enrichment <- m6A_enrichment %>% mutate(fdr.fisher.twosided = p.adjust(fisher.p.twosided, method = "fdr", n = length(fisher.p.twosided)))

                write.table(m6A_enrichment, file = paste0("~/Downstream_analysis/m6A_Genomic_enrichment/", tis, ".m6A.", dir, ".pval.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
	}
}

