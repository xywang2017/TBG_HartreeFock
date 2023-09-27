#!/bin/bash -l 

# for flag in {K,}
# do 
#     julia proj/runBMLL.jl $flag 07 1 8 1.20 strain
# done

# julia proj/runHF_server.jl 1 10 -3.1 ${2} 07 ${1} 1.20 strain symmetric

# julia proj/runHF_server.jl 1 8 -2.25 ${2} 07 ${1} 1.20 strain symmetric


julia proj/runBMLL.jl K 07 1 2 ${1} strain
julia proj/runBMLL.jl K 07 5 11 ${1} strain
julia proj/runBMLL.jl K 07 4 9 ${1} strain
julia proj/runBMLL.jl K 07 3 7 ${1} strain
julia proj/runBMLL.jl K 07 5 12 ${1} strain
# julia proj/runBMLL.jl K 07 2 5 ${1} strain
julia proj/runBMLL.jl K 07 3 8 ${1} strain
julia proj/runBMLL.jl K 07 4 11 ${1} strain
julia proj/runBMLL.jl K 07 1 3 ${1} strain
julia proj/runBMLL.jl K 07 3 10 ${1} strain
julia proj/runBMLL.jl K 07 2 7 ${1} strain
julia proj/runBMLL.jl K 07 3 11 ${1} strain
julia proj/runBMLL.jl K 07 1 4 ${1} strain
julia proj/runBMLL.jl K 07 2 9 ${1} strain
julia proj/runBMLL.jl K 07 1 5 ${1} strain
julia proj/runBMLL.jl K 07 2 11 ${1} strain
julia proj/runBMLL.jl K 07 1 6 ${1} strain
julia proj/runBMLL.jl K 07 1 7 ${1} strain
# julia proj/runBMLL.jl K 07 1 8 ${1} strain
julia proj/runBMLL.jl K 07 1 9 ${1} strain