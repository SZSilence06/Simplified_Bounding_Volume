#!/bin/bash

# usage : run file_path max_error sample_radius
# example : run bunny_2.obj 0.05 0.05
# to use this script, copy it to the binary directory and run it from that directory.

file_path=$1
max_error=$2
sample_radius=$3
file_name=$(echo $file_path | cut -f 1 -d '.')
result_dir="${file_name}_-e_${max_error}_-r_${sample_radius}"

cmd_1="./sbvgen --shell -s ${file_path} -e ${max_error} -r ${sample_radius}"
${cmd_1}
cmd_2="./sbvgen --inner ${result_dir}/inner_shell.vtk --outer ${result_dir}/outer_shell.vtk -d ${result_dir} -t" 
echo ${cmd_2}
${cmd_2}

