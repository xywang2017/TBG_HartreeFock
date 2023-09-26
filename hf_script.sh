#!/bin/bash -l 

for flag in {K,}
do 
    julia proj/runBMLL.jl $flag 07 3 7 1.20 strain
done

# julia proj/runHF_server.jl 2 5 -3.4 ${2} 07 ${1} 1.32 strain symmetric

# julia proj/runHF_server.jl 1 8 -1.375 ${2} 07 ${1} 1.32 strain symmetric