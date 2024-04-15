Phe1=$1
Phe2=$2

if [ $Phe1 = "m6A" ]; then
	Phe1_s="m6A-"
else
	Phe1_s=$Phe1
fi

if [ $Phe2 = "m6A" ]; then
        Phe2_s="m6A-"
else
        Phe2_s=$Phe2
fi

##generate .esi file for exposure modalities
for tis in Brain Lung Muscle Heart; do
	for i in `seq 22`; do
		awk -F '\t' '{if(NR==FNR){a[$0]=$0}else{if(a[$2]==$2){print $0}}}' \
			<(awk '{if(NR==FNR){a[$1]=$1}else{if(a[$1]==$1){print $2}}}' ~/2SMR/${Phe1}2${Phe2}/${tis}_${Phe1}2${Phe2}_pairs/${tis}.${Phe1}2${Phe2}_chr${i}.txt ~/Database/${tis}_${Phe1_s}QTL_signif/${tis}.${Phe1_s}QTL.signif_chr${i}.txt |sort -u) ~/Database/GTEx_v8/maf/chr${i}.frq |sed 's/:/\t/g' |\
			awk -v OFS='\t' '{print $1,$2":"$3":"$4":"$5,0,$3,$6,$7,$8}' > ~/SMR/${Phe1}2${Phe2}/${tis}_${Phe1_s}QTL_exp_besd/${tis}.${Phe1_s}QTL.exp_chr${i}.esi
	done
done

##generate .esi file for outcome modalities
for tis in Brain Lung Muscle Heart; do
	for i in `seq 22`; do
		awk -F '\t' '{if(NR==FNR){a[$0]=$0}else{if(a[$2]==$2){print $0}}}' \
			<(awk '{if(NR==FNR){a[$2]=$2}else{if(a[$1]==$1){print $2}}}' ~/2SMR/${Phe1}2${Phe2}/${tis}_${Phe1}2${Phe2}_pairs/${tis}.${Phe1}2${Phe2}_chr${i}.txt ~/Database/${tis}_${Phe2_s}QTL_all/${tis}.${Phe2_s}QTL.all_chr${i}.txt |sort -u) ~/Database/GTEx_v8/maf/chr${i}.frq |sed 's/:/\t/g' |\
			awk -v OFS='\t' '{print $1,$2":"$3":"$4":"$5,0,$3,$6,$7,$8}' > ~/SMR/${Phe1}2${Phe2}/${tis}_${Phe2_s}QTL_out_besd/${tis}.${Phe2_s}QTL.out_chr${i}.esi
	done
done

##generate .epi file for exposure modalities
for tis in Brain Lung Muscle Heart; do
	for i in `seq 22`; do
		awk '{if(NR==FNR){a[$1]=$1}else{if(a[$1]==$1){print $1}}}' ~/2SMR/${Phe1}2${Phe2}/${tis}_${Phe1}2${Phe2}_pairs/${tis}.${Phe1}2${Phe2}_chr${i}.txt ~/Database/${tis}_${Phe1_s}QTL_signif/${tis}.${Phe1_s}QTL.signif_chr${i}.txt |sort -u |sed 's/_/\t/g' |\
			awk -v OFS='\t' '{printf "%s\t%s\t%d\t%.0f\t%s\t%s\n",$1,$1"_"$2"_"$3,0,($2+$3)/2,$1"_"$2"_"$3,"+"}' > ~/SMR/${Phe1}2${Phe2}/${tis}_${Phe1_s}QTL_exp_besd/${tis}.${Phe1_s}QTL.exp_chr${i}.epi
	done
done

##generate .epi file for outcome modalities
for tis in Brain Lung Muscle Heart; do
	for i in `seq 22`; do
		awk '{if(NR==FNR){a[$2]=$2}else{if(a[$1]==$1){print $1}}}' ~/2SMR/${Phe1}2${Phe2}/${tis}_${Phe1}2${Phe2}_pairs/${tis}.${Phe1}2${Phe2}_chr${i}.txt ~/Database/${tis}_${Phe2_s}QTL_all/${tis}.${Phe2_s}QTL.all_chr${i}.txt |sort -u |sed 's/_/\t/g' |\
			awk -v OFS='\t' '{printf "%s\t%s\t%d\t%.0f\t%s\t%s\n",$1,$1"_"$2"_"$3,0,($2+$3)/2,$1"_"$2"_"$3,"+"}' > ~/SMR/${Phe1}2${Phe2}/${tis}_${Phe2_s}QTL_out_besd/${tis}.${Phe2_s}QTL.out_chr${i}.epi
	done
done

##generate .flist file for exposure modalities
for tis in Brain Lung Muscle Heart; do
	for i in `seq 22`; do
		awk -v tis=$tis -v Phe1=$Phe1 -v Phe2=$Phe2 -v Phe1_s=$Phe1_s -v OFS='\t' 'BEGIN{print "Chr","ProbeID","GeneticDistance","ProbeBp","Gene","Orientation","PathOfEsd"}{print $0,"~/SMR/"Phe1"2"Phe2"/"tis"_"Phe1_s"QTL_exp_esd/"$2".esd"}' ~/SMR/${Phe1}2${Phe2}/${tis}_${Phe1_s}QTL_exp_besd/${tis}.${Phe1_s}QTL.exp_chr${i}.epi > ~/SMR/${Phe1}2${Phe2}/${tis}_${Phe1_s}QTL_exp_besd/${tis}.${Phe1_s}QTL.exp_chr${i}.flist
	done
done

##generate .flist file for outcome modalities
for tis in Brain Lung Muscle Heart; do
	for i in `seq 22`; do
		awk -v tis=$tis -v Phe1=$Phe1 -v Phe2=$Phe2 -v Phe2_s=$Phe2_s -v OFS='\t' 'BEGIN{print "Chr","ProbeID","GeneticDistance","ProbeBp","Gene","Orientation","PathOfEsd"}{print $0,"~/SMR/"Phe1"2"Phe2"/"tis"_"Phe2_s"QTL_out_esd/"$2".esd"}' ~/SMR/${Phe1}2${Phe2}/${tis}_${Phe2_s}QTL_out_besd/${tis}.${Phe2_s}QTL.out_chr${i}.epi > ~/SMR/${Phe1}2${Phe2}/${tis}_${Phe2_s}QTL_out_besd/${tis}.${Phe2_s}QTL.out_chr${i}.flist
	done
done

##generate .besd file for exposure modalities
for tis in Brain Lung Muscle Heart; do
	for i in `seq 22`; do
		smr --eqtl-flist ~/SMR/${Phe1}2${Phe2}/${tis}_${Phe1_s}QTL_exp_besd/${tis}.${Phe1_s}QTL.exp_chr${i}.flist --make-besd-dense --out ~/SMR/${Phe1}2${Phe2}/${tis}_${Phe1_s}QTL_exp_besd/${tis}.${Phe1_s}QTL.exp_chr${i}
	done
done

##generate .besd file for outcome modalities
for tis in Brain Lung Muscle Heart; do
	for i in `seq 22`; do
		smr --eqtl-flist ~/SMR/${Phe1}2${Phe2}/${tis}_${Phe2_s}QTL_out_besd/${tis}.${Phe2_s}QTL.out_chr${i}.flist --make-besd-dense --out ~/SMR/${Phe1}2${Phe2}/${tis}_${Phe2_s}QTL_out_besd/${tis}.${Phe2_s}QTL.out_chr${i}
	done
done

##run smr cmd, Phe1 as exposure and Phe2 as outcome
for tis in Brain Lung Muscle Heart; do
	for i in `seq 22`; do
		smr --bfile ~/Database/GTEx_v8/bfile/chr${i} --beqtl-summary ~/SMR/${Phe1}2${Phe2}/${tis}_${Phe1_s}QTL_exp_besd/${tis}.${Phe1_s}QTL.exp_chr${i} --beqtl-summary ~/SMR/${Phe1}2${Phe2}/${tis}_${Phe2_s}QTL_out_besd/${tis}.${Phe2_s}QTL.out_chr${i} --peqtl-smr 0.99 --peqtl-heidi 0.99 --heidi-min-m 1 --out ~/SMR/${tis}_${Phe1}2${Phe2}_smr/chr${i} --smr-multi --ld-multi-snp 0.01
	done
done
