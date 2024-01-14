library(data.table)
library(dplyr)
library(coloc)
library(foreach)
library(doParallel)

coloc_abf <- function(i){
  Phe1_coloc_dat <- Phe1_dat %>% filter(Phenotype == coloc_pairs[i,exposure])
  Phe2_coloc_dat <- Phe2_dat %>% filter(Phenotype == coloc_pairs[i,outcome])
  coloc_dat <- inner_join(Phe1_coloc_dat, Phe2_coloc_dat, by = "SNP", suffix = c(paste0("_", Phe1), "_m6A"))
  coloc_res <- coloc.abf(dataset1 = list(snp = coloc_dat$SNP, beta = coloc_dat[[paste0("beta_", Phe1)]], varbeta = coloc_dat[[paste0("varbeta_", Phe1)]],
					 MAF = coloc_dat[[paste0("MAF_", Phe1)]], N = n[[Phe1]][tis], type = "quant"),
			 dataset2 = list(snp = coloc_dat$SNP, beta = coloc_dat$beta_m6A, varbeta = coloc_dat$varbeta_m6A,
					 MAF = coloc_dat$MAF_m6A, N = n[[Phe2]][tis], type = "quant"))
  coloc_results_file <- paste0("~/Coloc/", Phe1, "2m6A/", tis, "_", Phe1, "2m6A_results/", coloc_pairs[i,exposure], "-", coloc_pairs[i,outcome], ".coloc.results")
  coloc_results <- coloc_res$results %>% rename(SNP = snp) %>% arrange(desc(SNP.PP.H4))
  write.table(coloc_results, file = coloc_results_file, quote = FALSE, sep = "\t", row.names = FALSE)
  coloc_combine_res <- data.frame(exposure = coloc_pairs[i,exposure], outcome = coloc_pairs[i,outcome],
				  lead_SNP = coloc_results[1,"SNP"], t(as.data.frame(coloc_res$summary)))
  return(coloc_combine_res)
}

##Colocalization analyses for the same epigenome to m6A MR pairs as further validation
args <- commandArgs(T)
tis <- args[1]
n_chr <- args[2]
Phe1 <- args[3] ##"ha" or "me"

Phe2 <- "m6A"
##list of sample sizes for m6A-QTLs, haQTLs, and mQTLs in four tissues
n <- list("m6A"=c("Brain"=53, "Lung"=32, "Muscle"=40, "Heart"=40),
	  "ha"=c("Brain"=113, "Lung"=66, "Muscle"=108, "Heart"=100),
	  "me"=c("Lung"=190, "Muscle"=42))
##list of degrees of freedom for m6A-QTLs and haQTLs in four tissues, used for calculating se
df <- list("m6A"=c("Brain"=44, "Lung"=28, "Muscle"=33, "Heart"=33),
           "ha"=c("Brain"=96, "Lung"=62, "Muscle"=96, "Heart"=93))

chr_file <- paste0("~/2SMR/", Phe1, "2m6A/", tis, "_", Phe1, "2m6A_mr/chr", n_chr, ".mr.res")
Phe1_file <- paste0("~/Database/", tis, "_", Phe1, "QTL_all/", tis, ".", Phe1, "QTL.all_chr", n_chr, ".txt")
Phe2_file <- paste0("~/Database/", tis, "_m6A-QTL_all/", tis, ".m6A-QTL.all_chr", n_chr, ".txt")
maf_file <- paste0("~/Database/GTEx_v8/maf/chr", n_chr, ".frq")
out_file <- paste0("~/Coloc/", Phe1, "2m6A/", tis, "_", Phe1, "2m6A_coloc/chr", n_chr, ".coloc.res")

coloc_pairs <- fread(chr_file, sep = "\t", header = TRUE) %>% select(exposure, outcome)
maf <- fread(maf_file, sep = "\t", header = TRUE) %>% select(SNP, A1, A2, MAF)
if(Phe1 == "ha"){
	Phe1_dat <- fread(Phe1_file, sep = "\t", header = FALSE, col.names = c("Phenotype", "SNP", "beta", "effect_allele", "other_allele", "pvalues")) %>%
		inner_join(y = maf, by = "SNP") %>%
		mutate(varbeta = abs(beta/qt(pvalues/2, df[[Phe1]][tis]))*abs(beta/qt(pvalues/2, df[[Phe1]][tis])),
		       beta = ifelse(effect_allele == A1, beta, -beta))
}else{
	Phe1_dat <- fread(Phe1_file, sep = "\t", header = FALSE, col.names = c("Phenotype", "SNP", "beta", "se", "effect_allele", "other_allele", "MAF", "pvalues")) %>%
		inner_join(y = maf, by = "SNP") %>%
		mutate(varbeta = se*se, beta = ifelse(effect_allele == A1, beta, -beta))
}
Phe2_dat <- fread(Phe2_file, sep = "\t", header = FALSE, col.names = c("Phenotype", "SNP", "beta", "effect_allele", "other_allele", "pvalues")) %>%
	inner_join(y = maf, by = "SNP") %>%
	mutate(varbeta = abs(beta/qt(pvalues/2, df[[Phe2]][tis]))*abs(beta/qt(pvalues/2, df[[Phe2]][tis])),
	       beta = ifelse(effect_allele == A1, beta, -beta))

cl <- makeCluster(12)
registerDoParallel(cl)

res <- foreach(i = 1:nrow(coloc_pairs), .combine = rbind, .packages = c("dplyr", "coloc"), .errorhandling = "pass") %dopar% {
	coloc_abf(i)
}
write.table(res, file = out_file, quote = FALSE, sep = "\t", row.names = FALSE)

stopImplicitCluster()
stopCluster(cl)
