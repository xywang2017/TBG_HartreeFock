#!/bin/bash -l 

julia proj/runHF_server.jl 1 8 -3.125 random 07 ${1} 1.05 strain symmetric ${2} ${3} ${4}
