#!/bin/bash


#IN1=MoS2.sc.in
#OUT1=MoS2.sc.out

#IN2=MoS2.nsc.in
#OUT2=MoS2.nsc.out

#IN3=bands.in
#OUT3=bands.out


echo 'running job 1'
#mpirun -np 8 pw.x < $IN1 > $OUT1
mpirun -np 4 pw.x -i MoS2.sc.in -northo 8 > MoS2.sc.out

echo 'running job 2'
#mpirun -np 8 pw.x -northo 8  < $IN2 > $OUT2
mpirun -np 4 pw.x -i MoS2.nsc.in -northo 8  > MoS2.nsc.out

echo 'running job 3'
mpirun -np 4 bands.x -i bands.in -northo 8 > bands.out

echo 'done!'
