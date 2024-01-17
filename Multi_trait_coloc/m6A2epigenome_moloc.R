library(data.table)
library(dplyr)
library(moloc)
library(foreach)
library(doParallel)

moloc <- function(i){
	Phe1_moloc_dat <- Phe1_dat %>% filter(Phenotype == moloc_pairs[i, exposure])
	Phe2_moloc_dat <- Phe2_dat %>% filter(Phenotype == moloc_pairs[i, outcome])
	com_SNP <- Reduce(intersect, list(Phe1_moloc_dat$SNP, Phe2_moloc_dat$SNP, gwas_dat$SNP))
	if(length(com_SNP) > 1 && sort(gwas_dat[SNP %in% com_SNP,]$PVAL)[1] < 1e-4){
		listData <- list(Phe1_moloc_dat[SNP %in% com_SNP,], Phe2_moloc_dat[SNP %in% com_SNP,], gwas_dat[SNP %in% com_SNP,])
		moloc_res <- moloc_test(listData, prior_var = "default", save.SNP.info = FALSE)
		c_res <- c(t(moloc_res$best_snp))
		names(c_res) <- c("pp_a", "best_SNP_a", "pp_b", "best_SNP_b", "pp_ab", "best_SNP_ab", "pp_c", "best_SNP_c", "pp_ac", "best_SNP_ac", "pp_bc", "best_SNP_bc", "pp_abc", "best_SNP_abc")
		moloc_combine_res0 <- data.frame(exposure = moloc_pairs[i, exposure], outcome = moloc_pairs[i, outcome], nsnps = moloc_res$nsnps, t(c_res))
		return(moloc_combine_res0)
	}
}

##perform multi-trait colocalization analyses between significant m6A to epigenome MR pairs and GWAS
args <- commandArgs(T)
tis <- args[1]
nchr <- args[2]
Phe2 <- args[3] ##"ha" or "me"
gwasid <- args[4]

Phe1 <- "m6A"
##list of sample sizes for m6A-QTLs, haQTLs, and mQTLs in four tissues
n <- list("m6A"=c("Brain"=53, "Lung"=32, "Muscle"=40, "Heart"=40),
          "ha"=c("Brain"=113, "Lung"=66, "Muscle"=108, "Heart"=100),
          "me"=c("Lung"=190, "Muscle"=42))
##list of degrees of freedom for m6A-QTLs and haQTLs in four tissues, used for calculating se
df <- list("m6A"=c("Brain"=44, "Lung"=28, "Muscle"=33, "Heart"=33),
           "ha"=c("Brain"=96, "Lung"=62, "Muscle"=96, "Heart"=93))

chr_file <- paste0("~/Moloc/m6A2", Phe2, "/", tis, "_m6A2", Phe2, "_pairs/", tis, ".m6A2", Phe2, "_chr", nchr, ".txt")
Phe1_file <- paste0("~/Database/", tis, "_m6A-QTL_all/", tis, ".m6A-QTL.all_chr", nchr, ".txt")
Phe2_file <- paste0("~/Database/", tis, "_", Phe2, "QTL_all/", tis, ".", Phe2, "QTL.all_chr", nchr, ".txt")
maf_file <- paste0("~/Database/GTEx_v8/maf/chr", nchr, ".frq")
gwas_meta_file <- "~/Database/GWAS/GWAS_metadata.txt"
gwas_sum_file <- paste0("~/Moloc/sumstats/", gwasid, "/", gwasid, "_chr", nchr, ".sumstats")
out_file <- paste0("~/Moloc/m6A2", Phe2, "/", tis, "/", gwasid, "/chr", nchr, ".moloc.pp")

##read mol-QTLs data
moloc_pairs <- fread(chr_file, sep = "\t", header = FALSE) %>% dplyr::select(c(1,2,8,15)) %>% dplyr::rename(exposure = V1, outcome = V2, fdr = V8, PP4 = V15)
maf <- fread(maf_file, sep = "\t", header = TRUE) %>% dplyr::select(SNP, A1, A2, MAF)
Phe1_dat <- fread(Phe1_file, sep = "\t", header = FALSE, col.names = c("Phenotype", "SNP", "BETA", "effect_allele", "other_allele", "pval")) %>%
	inner_join(y = maf, by = "SNP") %>% mutate(SE = abs(BETA/qt(pval/2, df[[Phe1]][tis])), BETA = ifelse(effect_allele == A1, BETA, -BETA), N = n[[Phe1]][tis])
if(Phe2 == "ha"){
	Phe2_dat <- fread(Phe2_file, sep = "\t", header = FALSE, col.names = c("Phenotype", "SNP", "BETA", "effect_allele", "other_allele", "pval")) %>%
		inner_join(y = maf, by = "SNP") %>% mutate(SE = abs(BETA/qt(pval/2, df[[Phe2]][tis])), BETA = ifelse(effect_allele == A1, BETA, -BETA), N = n[[Phe2]][tis])
}else{
	Phe2_dat <- fread(Phe2_file, sep = "\t", header = FALSE, col.names = c("Phenotype", "SNP", "BETA", "SE", "effect_allele", "other_allele", "eaf", "pval")) %>%
		inner_join(y = maf, by = "SNP") %>% mutate(BETA = ifelse(effect_allele == A1, BETA, -BETA), N = n[[Phe2]][tis])
}

##processing GWAS summary data, harmonize effect allele and effect size
gwas_meta <- fread(gwas_meta_file, sep = "\t", header = FALSE, col.names = c("ID", "Ncases", "Nctrl", "N", "Type", "Description"))
gtex_file <- "~/Database/GTEx_v8/maf/chrall.txt"
gtex <- fread(gtex_file, sep = "\t", header = FALSE, col.names = c("CHR", "SNP", "Pos", "A1", "A2"))
gwas_dat <- fread(gwas_sum_file, sep = "\t", header = FALSE, col.names = c("CHR", "Pos", "A1_gwas", "A2_gwas", "Stats", "MAF", "PVAL")) %>% dplyr::select(-Stats) %>%
	inner_join(y = gtex, by = c("CHR", "Pos")) %>% mutate(A1_gwas = toupper(A1_gwas), A2_gwas = toupper(A2_gwas)) %>%
	filter((A1_gwas == A1 & A2_gwas == A2) | (A2_gwas == A1 & A1_gwas == A2)) %>% mutate(N = gwas_meta[ID == gwasid, N])

if(gwas_meta[ID == gwasid, Description] == "EAF"){
	gwas_dat <- gwas_dat %>% mutate(MAF = ifelse(A1_gwas == A1, MAF, 1 - MAF)) %>% select(SNP, A1, A2, MAF, PVAL, N)
}else{
	gwas_dat <- gwas_dat %>% select(SNP, A1, A2, MAF, PVAL, N)
}

if(gwas_meta[ID == gwasid, Type] == "binary"){
	gwas_dat <- gwas_dat %>% mutate(Ncases = gwas_meta[ID == gwasid, Ncases])
}

Phe1_dat <- Phe1_dat[complete.cases(Phe1_dat), ]
Phe2_dat <- Phe2_dat[complete.cases(Phe2_dat), ]
gwas_dat <- gwas_dat[complete.cases(gwas_dat), ]

cl <- makeCluster(12)
registerDoParallel(cl)

res <- foreach(i = 1:nrow(moloc_pairs), .combine = rbind, .packages = c("data.table", "dplyr", "moloc"), .errorhandling = "pass") %dopar% moloc(i)
write.table(res, file = out_file, quote = FALSE, sep = "\t", row.names = FALSE)
print(paste0(gwasid, " ", nchr, " moloc done"))

stopImplicitCluster()
stopCluster(cl)
