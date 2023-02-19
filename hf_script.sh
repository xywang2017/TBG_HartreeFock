#!/bin/bash -l 

for nu in {0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0}
do 
	julia proj/runHF_server.jl 1 4 $nu random 07 $1
done
