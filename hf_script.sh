# #!/bin/bash -l 

# for flag in {K,Kprime}
# do 
#     julia proj/runBMLL.jl $flag 07 1 4 1.20 strain & 
# done

julia proj/runHF_server.jl 1 10 -0.4 bm_cascade 07 1 1.20 strain
julia proj/runHF_server.jl 1 10 -1.3 bm_cascade 07 1 1.20 strain
julia proj/runHF_server.jl 1 10 -2.2 bm_cascade 07 1 1.20 strain
julia proj/runHF_server.jl 1 10 -3.1 bm_cascade 07 1 1.20 strain