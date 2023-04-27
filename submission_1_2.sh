#!/bin/bash
#SBATCH --job-name=work
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks=8
#SBATCH -A genacc_q
#SBATCH --mail-type="ALL" 
#SBATCH -t 72:00:00
#SBATCH --output=work%j.out
#SBATCH --error=work%j.err

module load julia

cd ~/server

julia proj/runHF_server.jl 1 2 0.000 random 07 1 &
julia proj/runHF_server.jl 1 2 0.500 random 07 1 &
julia proj/runHF_server.jl 1 2 1.000 random 07 1 &
julia proj/runHF_server.jl 1 2 1.500 random 07 1 &
julia proj/runHF_server.jl 1 2 2.000 random 07 1 &
julia proj/runHF_server.jl 1 2 2.500 random 07 1 &
julia proj/runHF_server.jl 1 2 3.000 random 07 1 &
julia proj/runHF_server.jl 1 2 3.500 random 07 1 &
wait
