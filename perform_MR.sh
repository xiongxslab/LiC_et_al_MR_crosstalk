perform_m6A2ha_MR(){
	for tis in Brain Lung Muscle Heart
	do
		for n in `seq 22`
		do
			Rscript m6A2epigenome_MR.R $tis $n ha
		done
	done
}

perform_m6A2me_MR(){
        for tis in Lung Muscle
        do
                for n in `seq 22`
                do
                        Rscript m6A2epigenome_MR.R $tis $n me
                done
        done
}

perform_ha2m6A_MR(){
        for tis in Brain Lung Muscle Heart
        do
                for n in `seq 22`
                do
                        Rscript epigenome2m6A_MR.R $tis $n ha
                done
        done
}

perform_me2m6A_MR(){
        for tis in Lung Muscle
        do
                for n in `seq 22`
                do
                        Rscript epigenome2m6A_MR.R $tis $n me
                done
        done
}

perform_m6A2ha_MR
#perform_m6A2me_MR
#perform_ha2m6A_MR
#perform_me2m6A_MR
