makesd_exp_m6A2ha(){
	for tis in Brain Lung Muscle Heart
	do
		for n in `seq 22`
		do
			Rscript makesd_exp.R $tis $n m6A ha
		done
	done
}

makesd_exp_m6A2me(){
        for tis in Lung Muscle
        do
                for n in `seq 22`
                do
                        Rscript makesd_exp.R $tis $n m6A me
                done
        done
}

makesd_exp_ha2m6A(){
        for tis in Brain Lung Muscle Heart
        do      
                for n in `seq 22`
                do      
                        Rscript makesd_exp.R $tis $n ha m6A
                done    
        done    
}

makesd_exp_me2m6A(){
        for tis in Lung Muscle
        do
                for n in `seq 22`
                do
                        Rscript makesd_exp.R $tis $n me m6A
                done
        done
}

makesd_out_m6A2ha(){
        for tis in Brain Lung Muscle Heart
        do      
                for n in `seq 22`
                do      
                        Rscript makesd_out.R $tis $n m6A ha
                done    
        done    
}

makesd_out_m6A2me(){
        for tis in Lung Muscle
        do
                for n in `seq 22`
                do
                        Rscript makesd_out.R $tis $n m6A me
                done
        done
}

makesd_out_ha2m6A(){
        for tis in Brain Lung Muscle Heart
        do
                for n in `seq 22`
                do
                        Rscript makesd_out.R $tis $n ha m6A
                done
        done
}

makesd_out_m6A2me(){
        for tis in Lung Muscle
        do
                for n in `seq 22`
                do
                        Rscript makesd_out.R $tis $n me m6A
                done
        done
}

makesd_exp_m6A2ha
#makesd_exp_m6A2me
#makesd_exp_ha2m6A
#makesd_exp_me2m6A
#makesd_out_m6A2ha
#makesd_out_m6A2me
#makesd_out_ha2m6A
#makesd_out_me2m6A
