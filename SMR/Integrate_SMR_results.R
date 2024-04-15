library(data.table)
library(dplyr)

for(tis in c("Brain", "Lung", "Muscle", "Heart")){
	for(dir in c("m6A2ha", "ha2m6A")){
		SMR_res <- data.frame()
		for(i in 1:22){
			SMR_res0 <- fread(paste0("~/SMR/", tis, "_", dir, "_smr/chr", i, ".msmr"), sep = "\t", header = TRUE)
			SMR_res <- bind_rows(SMR_res, SMR_res0)
		}

		SMR_res <- SMR_res %>% filter(!HEIDI_pval < 0.01)
		MR_res <- fread(paste0("~/2SMR/Results/", tis, ".", dir, ".mr.fdr.res"), sep = "\t", header = TRUE)
		Coloc_res <- fread(paste0("~/Coloc/Results/", tis, ".", dir, ".coloc.res"), sep = "\t", header = TRUE)
		Combine_res <- inner_join(x = MR_res, y = Coloc_res, by = c("exposure", "outcome")) %>%
			inner_join(y = SMR_res, by = c("exposure", "outcome"))

		write.table(Combine_res, file = paste0("~/SMR/Results/", tis, ".", dir, ".combine.res"), quote = FALSE, sep = "\t", row.names = FALSE)
	}
}

for(tis in c("Lung", "Muscle")){
	for(dir in c("m6A2me", "me2m6A")){
		SMR_res <- data.frame()
		for(i in 1:22){
			SMR_res0 <- fread(paste0("~/SMR/", tis, "_", dir, "_smr/chr", i, ".msmr"), sep = "\t", header = TRUE)
			SMR_res <- bind_rows(SMR_res, SMR_res0)
		}

		SMR_res <- SMR_res %>% filter(!HEIDI_pval < 0.01)
		MR_res <- fread(paste0("~/2SMR/Results/", tis, ".", dir, ".mr.fdr.res"), sep = "\t", header = TRUE)
		Coloc_res <- fread(paste0("~/Coloc/Results/", tis, ".", dir, ".coloc.res"), sep = "\t", header = TRUE)
		Combine_res <- inner_join(x = MR_res, y = Coloc_res, by = c("exposure", "outcome")) %>%
			inner_join(y = SMR_res, by = c("exposure", "outcome"))

		write.table(Combine_res, file = paste0("~/SMR/Results/", tis, ".", dir, ".combine.res"), quote = FALSE, sep = "\t", row.names = FALSE)
	}
}
