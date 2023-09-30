#!/bin/bash -l 

# for angle in {1.06,1.08,1.10,1.12,1.14,1.16,1.18,1.22,1.26,1.30,1.34,1.36}
# for angle in {1.05,1.20,1.24,1.28,1.32,1.38}
#     do 
#         for flag in {K,}
#         do 
#             for flag1 in {strain,nostrain}
#             do
#                 julia proj/runBMLL.jl $flag 07 2 5 $angle $flag1 
#             done
#         done
# done


for angle in {1.05,1.06,1.08,1.10,1.12,1.14,1.16,1.18,1.20,1.22,1.24,1.26,1.28,1.30,1.32,1.34,1.36,1.38}
    do 
        for flag1 in {strain,nostrain}
        do
            julia proj/runHF_server.jl 2 5 -3.4 bm_cascade 07 1 $angle $flag1 symmetric
        done
done
# julia proj/runHF_server.jl 3 11 -3.273 ${2} 07 ${1} 1.32 strain nosymmetric

# julia proj/runHF_server.jl 1 8 -2.25 ${2} 07 ${1} 1.20 strain symmetric


# julia proj/runBMLL.jl K 07 1 2 ${1} strain
# julia proj/runBMLL.jl K 07 5 11 ${1} strain
# julia proj/runBMLL.jl K 07 4 9 ${1} strain
# julia proj/runBMLL.jl K 07 3 7 ${1} strain
# julia proj/runBMLL.jl K 07 5 12 ${1} strain
# # julia proj/runBMLL.jl K 07 2 5 ${1} strain
# julia proj/runBMLL.jl K 07 3 8 ${1} strain
# julia proj/runBMLL.jl K 07 4 11 ${1} strain
# julia proj/runBMLL.jl K 07 1 3 ${1} strain
# julia proj/runBMLL.jl K 07 3 10 ${1} strain
# julia proj/runBMLL.jl K 07 2 7 ${1} strain
# julia proj/runBMLL.jl K 07 3 11 ${1} strain
# julia proj/runBMLL.jl K 07 1 4 ${1} strain
# julia proj/runBMLL.jl K 07 2 9 ${1} strain
# julia proj/runBMLL.jl K 07 1 5 ${1} strain
# julia proj/runBMLL.jl K 07 2 11 ${1} strain
# julia proj/runBMLL.jl K 07 1 6 ${1} strain
# julia proj/runBMLL.jl K 07 1 7 ${1} strain
# # julia proj/runBMLL.jl K 07 1 8 ${1} strain
# julia proj/runBMLL.jl K 07 1 9 ${1} strain