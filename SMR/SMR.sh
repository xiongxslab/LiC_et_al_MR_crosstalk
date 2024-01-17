##m6A to H3K27ac as example
##generate .esi file for m6A-QTLs as exposure
for tis in Brain Lung Muscle Heart
do 
	for i in `seq 22`
	do 
		awk -F '\t' '{if(NR==FNR){a[$0]=$0}else{if(a[$2]==$2){print $0}}}' \
			<(awk '{if(NR==FNR){a[$1]=$1}else{if(a[$1]==$1){print $2}}}' ~/2SMR/m6A2ha/${tis}_m6A2ha_pairs/${tis}.m6A2ha_chr${i}.txt ~/2SMR/m6A2ha/${tis}_m6A-QTL_signif/${tis}.m6A-QTL.signif_chr${i}.txt |sort -u) ~/Database/GTEx_v8/maf/chr${i}.frq |sed 's/:/\t/g' |\
			awk -v OFS='\t' '{print $1,$2":"$3":"$4":"$5,0,$3,$6,$7,$8}' > ~/SMR/m6A2ha/${tis}_m6A-QTL_exp_besd/${tis}.m6A-QTL.exp_chr${i}.esi
	done
done

##generate .epi file for m6A-QTLs as exposure
for tis in Brain Lung Muscle Heart
do
	for i in `seq 22`
	do
		awk '{if(NR==FNR){a[$1]=$1}else{if(a[$1]==$1){print $1}}}' ~/2SMR/m6A2ha/${tis}_m6A2ha_pairs/${tis}.m6A2ha_chr${i}.txt ~/2SMR/m6A2ha/${tis}_m6A-QTL_signif/${tis}.m6A-QTL.signif_chr${i}.txt |sort -u |sed 's/_/\t/g' |\
			awk -v OFS='\t' '{printf "%s\t%s\t%d\t%.0f\t%s\t%s\n",$1,$1"_"$2"_"$3,0,($2+$3)/2,$1"_"$2"_"$3,"+"}' |sed 's/^m6A:chr//g' > ~/SMR/m6A2ha/${tis}_m6A-QTL_exp_besd/${tis}.m6A-QTL.exp_chr${i}.epi
	done
done

##generate .flist file for m6A-QTLs as exposure
for tis in Brain Lung Muscle Heart
do
	for i in `seq 22`
	do
		awk -v tis=$tis -v OFS='\t' 'BEGIN{print "Chr","ProbeID","GeneticDistance","ProbeBp","Gene","Orientation","PathOfEsd"}{print $0,"~/SMR/m6A2ha/"tis"_m6A-QTL_exp_esd/"$2".esd"}' ~/SMR/m6A2ha/${tis}_m6A-QTL_exp_besd/${tis}.m6A-QTL.exp_chr${i}.epi > ~/SMR/m6A2ha/${tis}_m6A-QTL_exp_besd/${tis}.m6A-QTL.exp_chr${i}.flist
	done
done

##generate .besd file for m6A-QTLs as exposure
for tis in Brain Lung Muscle Heart
do
	for i in `seq 22`
	do
		smr --eqtl-flist ~/SMR/m6A2ha/${tis}_m6A-QTL_exp_besd/${tis}.m6A-QTL.exp_chr${i}.flist --make-besd-dense --out ~/SMR/m6A2ha/${tis}_m6A-QTL_exp_besd/${tis}.m6A-QTL.exp_chr${i}
	done
done

##run smr cmd, m6A-QTLs as exposure and haQTLs as outcome
for tis in Brain Lung Muscle Heart
do
	for i in `seq 22`
	do
		smr --bfile ~/Database/GTEx_v8/bfile/chr${i} --beqtl-summary ~/SMR/m6A2ha/${tis}_m6A-QTL_exp_besd/${tis}.m6A-QTL.exp_chr${i} --beqtl-summary ~/SMR/m6A2ha/${tis}_haQTL_out_besd/${tis}.haQTL.out_chr${i} --peqtl-smr 0.99 --peqtl-heidi 0.99 --heidi-min-m 1 --out ~/SMR/${tis}_m6A2ha_smr/chr${i} --smr-multi --ld-multi-snp 0.01
	done
done
