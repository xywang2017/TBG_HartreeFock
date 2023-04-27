#!/bin/bash
#SBATCH --job-name=work
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks=20
#SBATCH -A genacc_q
#SBATCH --mail-type="ALL" 
#SBATCH -t 72:00:00
#SBATCH --output=work%j.out
#SBATCH --error=work%j.err

module load julia

cd ~/server

julia proj/runHF_server.jl 1 8 0.000 random 07 1 &
julia proj/runHF_server.jl 1 8 0.125 random 07 1 &
julia proj/runHF_server.jl 1 8 0.250 random 07 1 &
julia proj/runHF_server.jl 1 8 0.375 random 07 1 &
julia proj/runHF_server.jl 1 8 0.500 random 07 1 &
julia proj/runHF_server.jl 1 8 1.000 random 07 1 &
julia proj/runHF_server.jl 1 8 1.125 random 07 1 &
julia proj/runHF_server.jl 1 8 1.250 random 07 1 &
julia proj/runHF_server.jl 1 8 1.375 random 07 1 &
julia proj/runHF_server.jl 1 8 1.500 random 07 1 &
julia proj/runHF_server.jl 1 8 2.000 random 07 1 &
julia proj/runHF_server.jl 1 8 2.125 random 07 1 &
julia proj/runHF_server.jl 1 8 2.250 random 07 1 &
julia proj/runHF_server.jl 1 8 2.375 random 07 1 &
julia proj/runHF_server.jl 1 8 2.500 random 07 1 &
julia proj/runHF_server.jl 1 8 3.000 random 07 1 &
julia proj/runHF_server.jl 1 8 3.125 random 07 1 &
julia proj/runHF_server.jl 1 8 3.250 random 07 1 &
julia proj/runHF_server.jl 1 8 3.375 random 07 1 &
julia proj/runHF_server.jl 1 8 3.500 random 07 1 &
wait
