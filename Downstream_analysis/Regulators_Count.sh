##Generate m6A-peak bed files
#all the m6A-peaks in the lung
cat ~/Database/Lung_m6A-QTL_all/* |cut -f 1 |sort -u |\
	sed 's/_/\t/g' |awk -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' |sort -k 1,1V -k 2,2n > ~/Downstream_analysis/Regulators/Lung.m6A_background.bed

#significant m6A-peaks in the lung
sed '1d' ~/2SMR/Results/Lung.m6A2me.mr.fdr.res |awk '{if($8<0.1){print $0}}' |cut -f 1 |sort -u |\
	sed 's/_/\t/g' |awk -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' |sort -k 1,1V -k 2,2n > ~/Downstream_analysis/Regulators/Lung.m6A_signif.m6A2me.bed

sed '1d' ~/2SMR/Results/Lung.me2m6A.mr.fdr.res |awk '{if($8<0.1){print $0}}' |cut -f 2 |sort -u |\
	sed 's/_/\t/g' |awk -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' |sort -k 1,1V -k 2,2n > ~/Downstream_analysis/Regulators/Lung.m6A_signif.me2m6A.bed

#all the DNAme sites in the lung, split by chromosomes
for i in `seq 22`; do
	cut -f 1 ~/Database/Lung_meQTL_all/Lung.meQTL.all_chr${i}.txt |sort -u |\
		awk -v OFS='\t' '{if(NR==FNR){a[$1]=$1;b[$1]=$2;c[$1]=$3;d[$1]=$4}else{if(a[$0]==$0){print b[$0],c[$0],d[$0],$0}}}' ~/Database/infinium-methylationepic.info - |\
		sort -k 1,1V -k 2,2n > ~/Downstream_analysis/Regulators/Lung_DNAme_split_bed/Lung.DNAme_background_chr${i}.bed
done

##Perform RBP and TF binding sites annotations for m6A-peaks and DNAme sites, respectively
Path=~/Downstream_analysis/Regulators
#RBP
for dir in m6A2me me2m6A; do
	signif_bed=${Path}/Lung.m6A_signif.${dir}.bed
	background_bed=${Path}/Lung.m6A_background.bed
	Anno_file=${Path}/human.bed.gz
	Out_file=${Path}/Lung.m6A_peak_RBP.${dir}.txt

	awk -v OFS='\t' '{if(NR==FNR){a[$1]=$1;b[$1]=$2}else{if(a[$1]==$1){print $1,b[$1],$2}else{print $1,0,$2}}}' \
		<(bedtools intersect -wo -a <(sed 's/^chr//g' ${signif_bed}) -b ${Anno_file} |cut -f 8 |sort |uniq -c |sed 's/^ *//g' |sed 's/ /\t/g' |awk '{print $2"\t"$1}') \
		<(bedtools intersect -wo -a <(sed 's/^chr//g' ${background_bed}) -b ${Anno_file} |cut -f 8 |sort |uniq -c |sed 's/^ *//g' |sed 's/ /\t/g' |awk '{print $2"\t"$1}') |sort -k 1,1V > ${Out_file}
done

#TF, calculate for each chromosome
for dir in m6A2me me2m6A; do
	for i in `seq 22`; do
		signif_bed=${Path}/Lung_DNAme_split_bed/DNAme_signif.${dir}_chr${i}.bed
		background_bed=${Path}/Lung_DNAme_split_bed/Lung.DNAme_background_chr${i}.bed
		Anno_file=${Path}/TFbs_split_bed/TF_binding_sites.chr${i}.bed.gz
		Out_file=${Path}/Lung_DNAme_split_bed/Lung.DNAme_TF.${dir}_chr${i}.txt
		
		awk -v OFS='\t' '{if(NR==FNR){a[$1]=$1;b[$1]=$2}else{if(a[$1]==$1){print $1,b[$1],$2}else{print $1,0,$2}}}' \
			<(bedtools intersect -wo -a ${signif_bed} -b ${Anno_file} |cut -f 8,9 |datamash -s -g 1 sum 2) \
			<(bedtools intersect -wo -a ${background_bed} -b ${Anno_file} |cut -f 8,9 |datamash -s -g 1 sum 2) |sort -k 1,1V > ${Out_file}
	done
done
