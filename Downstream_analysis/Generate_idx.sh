##Generate merged idx for m6A-peaks across tissues
cat ~/2SMR/Results/*.m6A2ha.mr.fdr.res |grep -v exposure |cut -f 1 |sort -u |sed 's/_/\t/g' |awk -v OFS='\t' '{print $1,$2,$3}' |sort -k 1,1V -k 2,2n |\
	bedtools merge -i - |awk '{print $0"\tm6A-"NR}' > ~/Downstream_analysis/Cross_tis_consistency/Merged.m6A_idx.m6A2ha.txt
cat ~/2SMR/Results/*.ha2m6A.mr.fdr.res |grep -v exposure |cut -f 2 |sort -u |sed 's/_/\t/g' |awk -v OFS='\t' '{print $1,$2,$3}' |sort -k 1,1V -k 2,2n |\
        bedtools merge -i - |awk '{print $0"\tm6A-"NR}' > ~/Downstream_analysis/Cross_tis_consistency/Merged.m6A_idx.ha2m6A.txt
cat ~/2SMR/Results/*.m6A2me.mr.fdr.res |grep -v exposure |cut -f 1 |sort -u |sed 's/_/\t/g' |awk -v OFS='\t' '{print $1,$2,$3}' |sort -k 1,1V -k 2,2n |\
        bedtools merge -i - |awk '{print $0"\tm6A-"NR}' > ~/Downstream_analysis/Cross_tis_consistency/Merged.m6A_idx.m6A2me.txt
cat ~/2SMR/Results/*.me2m6A.mr.fdr.res |grep -v exposure |cut -f 2 |sort -u |sed 's/_/\t/g' |awk -v OFS='\t' '{print $1,$2,$3}' |sort -k 1,1V -k 2,2n |\
        bedtools merge -i - |awk '{print $0"\tm6A-"NR}' > ~/Downstream_analysis/Cross_tis_consistency/Merged.m6A_idx.me2m6A.txt

##Generate merged idx for H3K27ac-peaks across tissues
cat ~/2SMR/Results/*.m6A2ha.mr.fdr.res |grep -v exposure |cut -f 2 |sort -u |sed 's/_/\t/g' |awk -v OFS='\t' '{print $1,$2,$3}' |sort -k 1,1V -k 2,2n |\
        bedtools merge -i - |awk '{print $0"\tha-"NR}' > ~/Downstream_analysis/Cross_tis_consistency/Merged.ha_idx.m6A2ha.txt
cat ~/2SMR/Results/*.ha2m6A.mr.fdr.res |grep -v exposure |cut -f 1 |sort -u |sed 's/_/\t/g' |awk -v OFS='\t' '{print $1,$2,$3}' |sort -k 1,1V -k 2,2n |\
        bedtools merge -i - |awk '{print $0"\tha-"NR}' > ~/Downstream_analysis/Cross_tis_consistency/Merged.ha_idx.ha2m6A.txt

##Index for m6A-peaks in each tissue
for tis in Brain Lung Muscle Heart; do
	sed '1d' ~/2SMR/Results/${tis}.m6A2ha.mr.fdr.res |cut -f 1 |sort -u |sed 's/_/\t/g' |awk -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' |sort -k 1,1V -k 2,2n |\
		bedtools intersect -wo -a - -b ~/Downstream_analysis/Cross_tis_consistency/Merged.m6A_idx.m6A2ha.txt |\
		awk '{print $8"\t"$4}' > ~/Downstream_analysis/Cross_tis_consistency/${tis}.m6A_idx.m6A2ha.txt
done

for tis in Brain Lung Muscle Heart; do
        sed '1d' ~/2SMR/Results/${tis}.ha2m6A.mr.fdr.res |cut -f 2 |sort -u |sed 's/_/\t/g' |awk -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' |sort -k 1,1V -k 2,2n |\
		bedtools intersect -wo -a - -b ~/Downstream_analysis/Cross_tis_consistency/Merged.m6A_idx.ha2m6A.txt |\
		awk '{print $8"\t"$4}' > ~/Downstream_analysis/Cross_tis_consistency/${tis}.m6A_idx.ha2m6A.txt
done

for tis in Lung Muscle; do
	sed '1d' ~/2SMR/Results/${tis}.m6A2me.mr.fdr.res |cut -f 1 |sort -u |sed 's/_/\t/g' |awk -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' |sort -k 1,1V -k 2,2n |\
		bedtools intersect -wo -a - -b ~/Downstream_analysis/Cross_tis_consistency/Merged.m6A_idx.m6A2me.txt |\
		awk '{print $8"\t"$4}' > ~/Downstream_analysis/Cross_tis_consistency/${tis}.m6A_idx.m6A2me.txt
done

for tis in Lung Muscle; do
	sed '1d' ~/2SMR/Results/${tis}.me2m6A.mr.fdr.res |cut -f 2 |sort -u |sed 's/_/\t/g' |awk -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' |sort -k 1,1V -k 2,2n |\
		bedtools intersect -wo -a - -b ~/Downstream_analysis/Cross_tis_consistency/Merged.m6A_idx.me2m6A.txt |\
		awk '{print $8"\t"$4}' > ~/Downstream_analysis/Cross_tis_consistency/${tis}.m6A_idx.me2m6A.txt
done

##Index for H3K27ac-peaks in each tissue
for tis in Brain Lung Muscle Heart; do
	sed '1d' ~/2SMR/Results/${tis}.m6A2ha.mr.fdr.res |cut -f 2 |sort -u |sed 's/_/\t/g' |awk -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' |sort -k 1,1V -k 2,2n |\
		bedtools intersect -wo -a - -b ~/Downstream_analysis/Cross_tis_consistency/Merged.ha_idx.m6A2ha.txt |\
		awk '{print $8"\t"$4}' > ~/Downstream_analysis/Cross_tis_consistency/${tis}.ha_idx.m6A2ha.txt
done

for tis in Brain Lung Muscle Heart; do
	sed '1d' ~/2SMR/Results/${tis}.ha2m6A.mr.fdr.res |cut -f 1 |sort -u |sed 's/_/\t/g' |awk -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' |sort -k 1,1V -k 2,2n |\
		bedtools intersect -wo -a - -b ~/Downstream_analysis/Cross_tis_consistency/Merged.ha_idx.ha2m6A.txt |\
		awk '{print $8"\t"$4}' > ~/Downstream_analysis/Cross_tis_consistency/${tis}.ha_idx.ha2m6A.txt
done
