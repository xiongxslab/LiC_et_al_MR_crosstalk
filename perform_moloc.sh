##Partition GWAS summary data into 22 chromosome files
map_file=~/Database/00-common_all.txt
cat ~/Database/GWAS/GWAS_metadata.txt |cut -f 1 |while read line
do
	mkdir -p ~/Moloc/sumstats/${line}
	zcat ~/Database/GWAS/Sumstats/${line}.sumstats.gz |\
		awk -v OFS='\t' '{if(NR==FNR){a[$3]=$3;b[$3]=$1;c[$3]=$2}else{if(a[$1]==$1){print b[$1],c[$1],$2,$3,$4,$5,$6}}}' ${map_file} - |\
		awk -v line=${line} -v tis=${tis} '{i=$1;print >> "~/Moloc/sumstats/"line/line"_chr"i".sumstats"}'
	
	echo "${line} partition done"
done
	
perform_m6A2ha_moloc(){
	cat ~/Database/GWAS/GWAS_metadata.txt |cut -f 1 |while read line
	do
        	for tis in Brain Lung Muscle Heart
        	do
                	for n in `seq 22`
                	do
                        	Rscript m6A2epigenome_moloc.R ${tis} ${n} ha ${line}
                	done
        	done
	done
}

perform_m6A2me_moloc(){
        cat ~/Database/GWAS/GWAS_metadata.txt |cut -f 1 |while read line
        do
                for tis in Lung Muscle
                do
                        for n in `seq 22`
                        do
                                Rscript m6A2epigenome_moloc.R ${tis} ${n} me ${line}
                        done
                done
        done    
}

perform_ha2m6A_moloc(){
        cat ~/Database/GWAS/GWAS_metadata.txt |cut -f 1 |while read line
        do
                for tis in Brain Lung Muscle Heart
                do
                        for n in `seq 22`
                        do
                                Rscript epigenome2m6A_moloc.R ${tis} ${n} ha ${line}
                        done
                done
        done
}

perform_me2m6A_moloc(){
        cat ~/Database/GWAS/GWAS_metadata.txt |cut -f 1 |while read line
        do
                for tis in Lung Muscle
                do
                        for n in `seq 22`
                        do
                                Rscript epigenome2m6A_moloc.R ${tis} ${n} me ${line}
                        done
                done
        done
}

perform_m6A2ha_moloc
#perform_m6A2me_moloc
#perform_ha2m6A_moloc
#perform_me2m6A_moloc
