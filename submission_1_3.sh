#!/bin/bash
#SBATCH --job-name=work
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks=12
#SBATCH -A genacc_q
#SBATCH --mail-type="ALL" 
#SBATCH -t 72:00:00
#SBATCH --output=work%j.out
#SBATCH --error=work%j.err

module load julia

cd ~/server

julia proj/runHF_server.jl 1 3 0.000 random 07 1 &
julia proj/runHF_server.jl 1 3 0.333 random 07 1 &
julia proj/runHF_server.jl 1 3 0.667 random 07 1 &
julia proj/runHF_server.jl 1 3 1.000 random 07 1 &
julia proj/runHF_server.jl 1 3 1.333 random 07 1 &
julia proj/runHF_server.jl 1 3 1.667 random 07 1 &
julia proj/runHF_server.jl 1 3 2.000 random 07 1 &
julia proj/runHF_server.jl 1 3 2.333 random 07 1 &
julia proj/runHF_server.jl 1 3 2.667 random 07 1 &
julia proj/runHF_server.jl 1 3 3.000 random 07 1 &
julia proj/runHF_server.jl 1 3 3.333 random 07 1 &
julia proj/runHF_server.jl 1 3 3.667 random 07 1 &
wait
