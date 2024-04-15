library(data.table)
library(dplyr)

##Integrate m6A-H3K27ac MR results for the four tissues
for(tis in c("Brain", "Lung", "Muscle", "Heart")){
	for(dir in c("m6A2ha", "ha2m6A")){
		MR_res <- data.frame()
		for(i in 1:22){
			MR_res0 <- fread(paste0("~/2SMR/", dir, "/", tis, "_", dir, "_mr/chr", i, ".mr.res"), sep = "\t", header = TRUE)
			MR_res <- bind_rows(MR_res, MR_res0)
		}
		
		MR_res <- MR_res %>% mutate(fdr = p.adjust(pval, method = "fdr", n = length(pval)))
		write.table(MR_res, file = paste0("~/2SMR/Results/", tis, ".", dir, ".mr.fdr.res"), quote = FALSE, sep = "\t", row.names = FALSE)
	}
}

##Integrate m6A-DNAme MR results for lung and muscle
for(tis in c("Lung", "Muscle")){
	for(dir in c("m6A2me", "me2m6A")){
		MR_res <- data.frame()
		for(i in 1:22){
			MR_res0 <- fread(paste0("~/2SMR/", dir, "/", tis, "_", dir, "_mr/chr", i, ".mr.res"), sep = "\t", header = TRUE)
			MR_res <- bind_rows(MR_res, MR_res0)
		}

		MR_res <- MR_res %>% mutate(fdr = p.adjust(pval, method = "fdr", n = length(pval)))
		write.table(MR_res, file = paste0("~/2SMR/Results/", tis, ".", dir, ".mr.fdr.res"), quote = FALSE, sep = "\t", row.names = FALSE)
	}
}
