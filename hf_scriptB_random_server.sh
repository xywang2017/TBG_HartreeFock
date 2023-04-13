#!/bin/bash -l 

# for nu in {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5}
# do 
# 	julia proj/runHF_server.jl 1 2 $nu random 07 $1 &
# done

for nu in {0.0,0.429,0.571,0.857,1.0,1.286,1.429,1.857,2.0,2.429,2.857,3.0,3.429}
do 
	julia proj/runHF_server.jl 3 7 $nu random 07 $1 &
done


for nu in {0.0,0.4,0.6,0.8,1.0,1.2,1.4,1.8,2.0,2.4,2.8,3.0,3.4}
do 
	julia proj/runHF_server.jl 2 5 $nu random 07 $1 &
done

for nu in {0.0,0.286,0.571,0.714,0.857,1.0,1.286,1.571,2.0,2.286,2.571,3.0,3.286}
do 
	julia proj/runHF_server.jl 2 7 $nu random 07 $1 &
done


# for nu in {0.0,0.333,0.667,1.0,1.333,1.667,2.0,2.333,2.667,3.0,3.333}
# do 
# 	julia proj/runHF_server.jl 1 3 $nu random 07 $1 &
# done

# for nu in {0.0,0.25,0.5,0.75,1.0,1.25,1.5,2.0,2.25,2.5,3.0,3.25}
# do 
# 	julia proj/runHF_server.jl 1 4 $nu random 07 $1 &
# done

# for nu in {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,2.0,2.2,2.4,3.0,3.2}
# do 
# 	julia proj/runHF_server.jl 1 5 $nu random 07 $1 &
# done

# for nu in {0.0,0.167,0.333,0.5,0.833,1.0,1.167,1.333,2.0,2.167,2.333,3.0,3.167}
# do 
# 	julia proj/runHF_server.jl 1 6 $nu random 07 $1 &
# done

# for nu in {0.0,0.125,0.25,0.375,0.875,1.0,1.125,1.25,2.0,2.125,2.25,3.0,3.125}
# do 
# 	julia proj/runHF_server.jl 1 8 $nu random 07 $1 &
# done