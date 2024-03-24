##Annotate m6A loci in all tested m6A-DNAme MR pairs using Homer
for tis in Lung Muscle; do
	cat <(sed '1d' ~/2SMR/Results/${tis}.m6A2me.mr.fdr.res |cut -f 1) <(sed '1d' ~/2SMR/Results/${tis}.me2m6A.mr.fdr.res |cut -f 2) \|
		sort -u |sed 's/_/\t/g' |awk -v OFS='\t' '{print $1"_"$2"_"$3,$1,$2,$3,"+"}' |sort -k 2,2V -k 3,3n > ~/Downstream_analysis/m6A_Genomic_enrichment/${tis}.m6A_peak.merged_tested.bed
	
	annotatePeaks.pl ~/Downstream_analysis/m6A_Genomic_enrichment/${tis}.m6A_peak.merged_tested.bed hg38 > ~/Downstream_analysis/m6A_Genomic_enrichment/${tis}.m6A_peak.merged_tested.anno.txt
done
