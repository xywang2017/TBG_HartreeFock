#!/bin/bash
#SBATCH --job-name=work
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks=18
#SBATCH -A genacc_q
#SBATCH --mail-type="ALL" 
#SBATCH -t 72:00:00
#SBATCH --output=work%j.out
#SBATCH --error=work%j.err

module load julia

cd ~/server

julia proj/runHF_server.jl 3 7 0.000 random 07 1 &
julia proj/runHF_server.jl 3 7 0.429 random 07 1 &
julia proj/runHF_server.jl 3 7 0.857 random 07 1 &
julia proj/runHF_server.jl 3 7 1.000 random 07 1 &
julia proj/runHF_server.jl 3 7 1.286 random 07 1 &
julia proj/runHF_server.jl 3 7 1.429 random 07 1 &
julia proj/runHF_server.jl 3 7 1.714 random 07 1 &
julia proj/runHF_server.jl 3 7 1.857 random 07 1 &
julia proj/runHF_server.jl 3 7 2.000 random 07 1 &
julia proj/runHF_server.jl 3 7 2.286 random 07 1 &
julia proj/runHF_server.jl 3 7 2.429 random 07 1 &
julia proj/runHF_server.jl 3 7 2.714 random 07 1 &
julia proj/runHF_server.jl 3 7 2.857 random 07 1 &
julia proj/runHF_server.jl 3 7 3.000 random 07 1 &
julia proj/runHF_server.jl 3 7 3.286 random 07 1 &
julia proj/runHF_server.jl 3 7 3.429 random 07 1 &
julia proj/runHF_server.jl 3 7 3.714 random 07 1 &
julia proj/runHF_server.jl 3 7 3.857 random 07 1 &
wait
