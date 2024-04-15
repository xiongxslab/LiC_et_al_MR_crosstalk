library(data.table)
library(dplyr)

for(tis in c("Brain", "Lung", "Muscle", "Heart")){
	for(dir in c("m6A2ha", "ha2m6A")){
		Coloc_res <- data.frame()
		for(i in 1:22){
			Coloc_res0 <- fread(paste0("~/Coloc/", dir, "/", tis, "_", dir, "_coloc/chr", i, ".coloc.res"), sep = "\t", header = TRUE)
			Coloc_res <- bind_rows(Coloc_res, Coloc_res0)
		}

		write.table(Coloc_res, file = paste0("~/Coloc/Results/", tis, ".", dir, ".coloc.res"), quote = FALSE, sep = "\t", row.names = FALSE)
	}
}

for(tis in c("Lung", "Muscle")){
	for(dir in c("m6A2me", "me2m6A")){
		Coloc_res <- data.frame()
		for(i in 1:22){
			Coloc_res0 <- fread(paste0("~/Coloc/", dir, "/", tis, "_", dir, "_coloc/chr", i, ".coloc.res"), sep = "\t", header = TRUE)
			Coloc_res <- bind_rows(Coloc_res, Coloc_res0)
		}

		write.table(Coloc_res, file = paste0("~/Coloc/Results/", tis, ".", dir, ".coloc.res"), quote = FALSE, sep = "\t", row.names = FALSE)
	}
}
