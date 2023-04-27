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

julia proj/runHF_server.jl 2 5 0.000 random 07 1 &
julia proj/runHF_server.jl 2 5 0.400 random 07 1 &
julia proj/runHF_server.jl 2 5 0.800 random 07 1 &
julia proj/runHF_server.jl 2 5 1.000 random 07 1 &
julia proj/runHF_server.jl 2 5 1.200 random 07 1 &
julia proj/runHF_server.jl 2 5 1.400 random 07 1 &
julia proj/runHF_server.jl 2 5 1.600 random 07 1 &
julia proj/runHF_server.jl 2 5 1.800 random 07 1 &
julia proj/runHF_server.jl 2 5 2.000 random 07 1 &
julia proj/runHF_server.jl 2 5 2.200 random 07 1 &
julia proj/runHF_server.jl 2 5 2.400 random 07 1 &
julia proj/runHF_server.jl 2 5 2.600 random 07 1 &
julia proj/runHF_server.jl 2 5 2.800 random 07 1 &
julia proj/runHF_server.jl 2 5 3.000 random 07 1 &
julia proj/runHF_server.jl 2 5 3.200 random 07 1 &
julia proj/runHF_server.jl 2 5 3.400 random 07 1 &
julia proj/runHF_server.jl 2 5 3.600 random 07 1 &
julia proj/runHF_server.jl 2 5 3.800 random 07 1 &
wait
