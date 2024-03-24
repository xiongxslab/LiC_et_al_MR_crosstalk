##Generate DNAme bed files for ChromHMM annotations
#all tested DNAme
for tis in Lung Muscle; do
	sed '1d' ~/2SMR/Results/${tis}.m6A2me.mr.fdr.res |cut -f 2 |sort -u |\
		awk -v OFS='\t' '{if(NR==FNR){a[$1]=$1;b[$1]=$2;c[$1]=$3;d[$1]=$4}else{if(a[$0]==$0){print b[$0],c[$0],d[$0],$0}}}' ~/Database/infinium-methylationepic.info - |\
		sort -k 1,1V -k 2,2n > ~/Downstream_analysis/Chrom_states_enrichment/${tis}.DNAme_tested.m6A2me.bed
done

for tis in Lung Muscle; do
	sed '1d' ~/2SMR/Results/${tis}.me2m6A.mr.fdr.res |cut -f 1 |sort -u |\
		awk -v OFS='\t' '{if(NR==FNR){a[$1]=$1;b[$1]=$2;c[$1]=$3;d[$1]=$4}else{if(a[$0]==$0){print b[$0],c[$0],d[$0],$0}}}' ~/Database/infinium-methylationepic.info - |\
		sort -k 1,1V -k 2,2n > ~/Downstream_analysis/Chrom_states_enrichment/${tis}.DNAme_tested.me2m6A.bed
done

#significant DNAme
for tis in Lung Muscle; do
	sed '1d' ~/2SMR/Results/${tis}.m6A2me.mr.fdr.res |awk '{if($8<0.1){print $0}}' |cut -f 2 |sort -u |\
		awk -v OFS='\t' '{if(NR==FNR){a[$1]=$1;b[$1]=$2;c[$1]=$3;d[$1]=$4}else{if(a[$0]==$0){print b[$0],c[$0],d[$0],$0}}}' ~/Database/infinium-methylationepic.info - |\
		sort -k 1,1V -k 2,2n > ~/Downstream_analysis/Chrom_states_enrichment/${tis}.DNAme_signif.m6A2me.bed
done

for tis in Lung Muscle; do
	sed '1d' ~/2SMR/Results/${tis}.me2m6A.mr.fdr.res |awk '{if($8<0.1){print $0}}' |cut -f 1 |sort -u |\
		awk -v OFS='\t' '{if(NR==FNR){a[$1]=$1;b[$1]=$2;c[$1]=$3;d[$1]=$4}else{if(a[$0]==$0){print b[$0],c[$0],d[$0],$0}}}' ~/Database/infinium-methylationepic.info - |\
		sort -k 1,1V -k 2,2n > ~/Downstream_analysis/Chrom_states_enrichment/${tis}.DNAme_signif.me2m6A.bed
done

##Perform ChromHMM annotations for significant and all tested DNAme sites
Path=~/Downstream_analysis/Chrom_states_enrichment
#Lung, use epigenome E096 as annotation
for dir in m6A2me me2m6A; do
	signif_bed=${Path}/Lung.DNAme_signif.${dir}.bed
	tested_bed=${Path}/Lung.DNAme_tested.${dir}.bed
	Anno_file=${Path}/E096_18_core_K27ac_hg38lift_mnemonics.bed.gz
	Out_file=${Path}/Lung.DNAme.${dir}.txt

	awk -v OFS='\t' '{if(NR==FNR){a[$1]=$1;b[$1]=$2}else{if(a[$1]==$1){print $1,b[$1],$2}else{print $1,0,$2}}}' \
		<(bedtools intersect -wo -a ${signif_bed} -b ${Anno_file} |cut -f 4,8 |sort -u |cut -f 2 |sort |uniq -c |sed 's/^ *//g' |sed 's/ /\t/g' |awk '{print $2"\t"$1}') \
		<(bedtools intersect -wo -a ${tested_bed} -b ${Anno_file} |cut -f 4,8 |sort -u |cut -f 2 |sort |uniq -c |sed 's/^ *//g' |sed 's/ /\t/g' |awk '{print $2"\t"$1}') |sort -k 1,1V > ${Out_file}
done

#Muscle, use epigenome E108 as annotation
for dir in m6A2me me2m6A; do
	signif_bed=${Path}/Muscle.DNAme_signif.${dir}.bed
	tested_bed=${Path}/Muscle.DNAme_tested.${dir}.bed
	Anno_file=${Path}/E108_18_core_K27ac_hg38lift_mnemonics.bed.gz
	Out_file=${Path}/Muscle.DNAme.${dir}.txt

	awk -v OFS='\t' '{if(NR==FNR){a[$1]=$1;b[$1]=$2}else{if(a[$1]==$1){print $1,b[$1],$2}else{print $1,0,$2}}}' \
		<(bedtools intersect -wo -a ${signif_bed} -b ${Anno_file} |cut -f 4,8 |sort -u |cut -f 2 |sort |uniq -c |sed 's/^ *//g' |sed 's/ /\t/g' |awk '{print $2"\t"$1}') \
		<(bedtools intersect -wo -a ${tested_bed} -b ${Anno_file} |cut -f 4,8 |sort -u |cut -f 2 |sort |uniq -c |sed 's/^ *//g' |sed 's/ /\t/g' |awk '{print $2"\t"$1}') |sort -k 1,1V > ${Out_file}
done
