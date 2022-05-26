LANG=en_US
let sample=1
for lr in $(seq 0.01 0.01 0.30)
do
	for deltachi in $(seq 1 5 100)
	do
		echo IIIIIIIIIIIIIIIIIIIIIIIII ${sample} IIIIIIIIIIIIIIIIIIIIIII
		mkdir out_${lr}_${deltachi}
		let sample++

		root -b -q 'GenEv_refit.C('$lr','$deltachi')'
		root -b -q dohist.C
	
		mv ${lr}_${deltachi} /home/piotrek/RefitNew/Results/out_${lr}_${deltachi}
		mv *.png /home/piotrek/RefitNew/Results/out_${lr}_${deltachi}
		mv Output_GenEv_1.root hist.root Refit.root /home/piotrek/RefitNew/Results/out_${lr}_${deltachi}
	done
done