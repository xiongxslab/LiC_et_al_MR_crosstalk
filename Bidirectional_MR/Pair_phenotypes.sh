for tis in Brain Lung Muscle Heart; do
	for n in `seq 22`; do
		python Pair_phenotypes.py $tis m6A ha $n
		python Pair_phenotypes.py $tis ha m6A $n
	done
done

for tis in Lung Muscle; do
	for n in `seq 22`; do
		python Pair_phenotypes.py $tis m6A me $n
		python Pair_phenotypes.py $tis me m6A $n
	done
done
