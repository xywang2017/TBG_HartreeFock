#!/bin/bash -l 

for seed in {3..5}
do
 	julia proj/runHF_server.jl 1 4 1.0 random 07 $seed &
done
# for q in {2,3,4,5,6,8,10}
# do 
# 	julia proj/runHF_server.jl 1 $q 0 flavor 07 $1
# done

# for nu in {0.0,0.5,1.0,1.5,2.0,2.5,3.0}
# do 
# 	julia proj/runHF_server.jl 1 2 $nu random 07 $1
# done

# for nu in {0.0,0.333,0.667,1.0,1.333,1.667,2.0,2.333,2.667,3.0}
# do 
# 	julia proj/runHF_server.jl 1 3 $nu random 07 $1
# done

# for nu in {0.0,0.25,0.5,0.75,1.0,2.0,2.25,2.5}
# do 
# 	julia proj/runHF_server.jl 1 4 $nu random 07 $1
# done

#for nu in {0.0,0.8,1.6,2.4,3.2}
#do 
#	julia proj/runHF_server.jl 1 5 $nu random 07 $1
#done


#for nu in {0.0,0.667,1.5,2.333,3.167}
#do 
#	julia proj/runHF_server.jl 1 6 $nu random 07 $1
#done

# for nu in {0.0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1.0,2.0,2.125,2.25,2.375,2.5}
# for nu in {1.125,1.25,1.375,1.5,1.625,1.75,1.875}
# for nu in {1.0,}
# do 
# 	julia --threads=1 proj/runHF_server.jl 1 8 $nu random 07 $1
# done

#for nu in {0.0,0.4,1.3,2.2,3.1}
#do 
#	julia proj/runHF_server.jl 1 10 $nu random 07 $1
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