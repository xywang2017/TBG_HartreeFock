#!/bin/bash -l 

# for nu in {0.0,0.5,1.0,1.75,2.0,2.25,2.5}
for nu in {0.25,}
do 
	julia proj/runHF_server.jl 1 4 $nu random 07 $1
done

#for nu in {0.0,0.8,1.6,2.4,3.2}
#do 
#	julia proj/runHF_server.jl 1 5 $nu random 07 $1
#done


#for nu in {0.0,0.667,1.5,2.333,3.167}
#do 
#	julia proj/runHF_server.jl 1 6 $nu random 07 $1
#done


#for nu in {0.0,0.4,1.3,2.2,3.1}
#do 
#	julia proj/runHF_server.jl 1 10 $nu random 07 $1
#done

#for nu in {0.0,0.5,1.375,2.25,3.125}
#do 
#	julia proj/runHF_server.jl 1 8 $nu random 07 $1
#done

# julia proj/runHF_server.jl 1 2 2.5 random 07 $1
# julia proj/runHF_server.jl 5 12 2.25 random 07 $1

# julia proj/runHF_server.jl 2 5 2.2 random 07 $1

# julia proj/runHF_server.jl 1 3 2.0 random 07 $1
# julia proj/runHF_server.jl 2 7 1.857 random 07 $1

# julia proj/runHF_server.jl 1 4 1.75 random 07 $1

# julia proj/runHF_server.jl 1 5 1.6 random 07 $1

# julia proj/runHF_server.jl 1 6 1.5 random 07 $1

# julia proj/runHF_server.jl 1 8 1.375 random 07 $1

# julia proj/runHF_server.jl 1 10 1.3 random 07 $1

# julia proj/runHF_server.jl 1 12 1.25 random 07 $1

# julia proj/runHF_server.jl 1 14 1.214 random 07 $1