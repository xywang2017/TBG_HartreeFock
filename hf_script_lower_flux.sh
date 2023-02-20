#!/bin/bash -l 

for nu in {0.0,0.8,0.6,2.4,3.2}
do 
	julia proj/runHF_server.jl 1 5 $nu random 07 $1
done

for nu in {0.0,0.667,1.5,2.333,3.167}
do 
	julia proj/runHF_server.jl 1 6 $nu random 07 $1
done

for nu in {0.0,0.5,1.375,2.25,3.125}
do 
	julia proj/runHF_server.jl 1 6 $nu random 07 $1
done