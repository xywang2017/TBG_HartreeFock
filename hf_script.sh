# #!/bin/bash -l 

# for flag in {K,Kprime}
# do 
#     julia proj/runBMLL.jl $flag 07 1 4 1.20 strain & 
# done

julia proj/runHF_server.jl 1 4 -1.75 bm_cascade 07 1 1.20 strain