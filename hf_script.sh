# #!/bin/bash -l 

# for flag in {K,}
# do 
#     julia proj/runBMLL.jl $flag 07 2 5 1.20 strain & 
# done

julia proj/runHF_server.jl 1 10 -2.2 bm_cascade 07 2 1.20 strain