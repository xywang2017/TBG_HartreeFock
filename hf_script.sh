# #!/bin/bash -l 

# for flag in {K,Kprime}
# do 
#     julia proj/runBMLL.jl $flag 07 1 5 1.20 strain & 
# done

julia proj/runHF_server.jl 1 5 -2.4 bm_cascade 07 1 1.20 strain