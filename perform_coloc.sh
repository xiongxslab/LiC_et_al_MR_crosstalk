perform_m6A2ha_coloc(){
	for tis in Brain Lung Muscle Heart
	do
		for n in `seq 22`
		do
			Rscript m6A2epigenome_coloc.R $tis $n ha
		done
	done
}

perform_m6A2me_coloc(){
        for tis in Lung Muscle
        do
                for n in `seq 22`
                do
                        Rscript m6A2epigenome_coloc.R $tis $n me
                done
        done
}

perform_ha2m6A_coloc(){
        for tis in Brain Lung Muscle Heart
        do
                for n in `seq 22`
                do
                        Rscript epigenome2m6A_coloc.R $tis $n ha
                done
        done
}

perform_me2m6A_coloc(){
        for tis in Lung Muscle
        do
                for n in `seq 22`
                do
                        Rscript epigenome2m6A_coloc.R $tis $n me
                done
        done
}

perform_m6A2ha_coloc
#perform_m6A2me_coloc
#perform_ha2m6A_coloc
#perform_me2m6A_coloc
