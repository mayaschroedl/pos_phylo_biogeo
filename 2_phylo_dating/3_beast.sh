#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 4G
#SBATCH -c 1
#SBATCH --time 2-0

../../beast/bin/beast most_diff_concat.xml