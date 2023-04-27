using Printf 

ϕs = [1//8;1//6;1//5;1//4;2//7;1//3;2//5;3//7;1//2]
for ϕ in ϕs 
    p,q = numerator(ϕ), denominator(ϕ)
    νs=unique(sort(round.([s+t*p/q for s in 0:3 for t in 0:4],digits=3)))
    νs = νs[νs .< 4]
    num = length(νs)
    open("submission_$(p)_$(q).sh","w") do file 
        write(file,"#!/bin/bash\n")
        write(file,"#SBATCH --job-name=work\n")
        write(file,"#SBATCH --ntasks-per-core=1\n")
        write(file,"#SBATCH --ntasks=$(num)\n")
        write(file,"#SBATCH -A genacc_q\n")
        write(file,"#SBATCH --mail-type=\"ALL\" \n")
        write(file,"#SBATCH -t 72:00:00\n")
        write(file,"#SBATCH --output=work%j.out\n")
        write(file,"#SBATCH --error=work%j.err\n\n")
        write(file,"module load julia\n\n")
        write(file,"cd ~/server\n\n")
        for i in eachindex(νs)
            write(file,@sprintf "julia proj/runHF_server.jl %d %d %.3f %s %s %d &\n" p q νs[i] "random" "07" 1)
        end
        write(file,"wait\n")
    end
end