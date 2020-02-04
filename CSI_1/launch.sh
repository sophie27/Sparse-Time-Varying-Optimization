#!/bin/bash
export LC_ALL="C";
T_zero="$(($(date +%s)))"
N_RUNS=$(seq 1 200)
N_PARALLEL_TASKS=5
N=20
M=12
K=2
#np=0.1 #SNR=20dB if u= Standard Gaussian
np=0.0562 #SNR=25dB if u= Standard Gaussian
#np=0

rm csi*.dat

RUN_SESSION=0

for RUN in $N_RUNS; do
	if [ $RUN_SESSION -eq $N_PARALLEL_TASKS ]; then	
		while [ 1 -eq 1 ]; do
			sleep 1
			RUNNING_PROCESSES=$(ps aux | grep dynreg | grep -v grep | wc -l)
			if [ $RUNNING_PROCESSES -lt $N_PARALLEL_TASKS ]; then
				RUN_SESSION=$(($RUNNING_PROCESSES-1))
				break
			fi
		done
	fi
	./dynreg $N $K $M $RUN $np > "csi${RUN}.dat" &
	RUN_SESSION=$(($RUN_SESSION + 1))
done

while [ 1 -eq 1 ]; do
		sleep 1
		RUNNING_PROCESSES=$(ps aux | grep dynreg | grep -v grep | wc -l)
			if [ $RUNNING_PROCESSES -eq 0 ]; then
				break
			fi
done

rm csi.dat
nohup matlab -nodisplay < merge.m > out.out
gnuplot plot.gnu


T="$(($(date +%s)-T_zero))"
echo $T > Time_Elapsed.txt
#

