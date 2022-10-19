#!/bin/bash
#ref:https://github.com/Teichlab/cellphonedb
input_road=$1
meta_file=$2
count_file=$3
#subtype_ratio=$4
project_name=$4
output_road=$5

cellphonedb=/mnt/data/chenwei/software/miniconda2_new/envs/cpdb/bin/cellphonedb
$cellphonedb method statistical_analysis \
              $input_road/${meta_file}  $input_road/${count_file} \
              --counts-data=gene_name \
              --iterations=1000 \
              --threads=30 \
              --project-name=$project_name\
              --output-path=$output_road
mkdir $output_road/$project_name
mv ./out/*  $output_road/$project_name

$cellphonedb plot dot_plot \
            --means-path=${output_road}/${project_name}/means.txt \
            --pvalues-path=${output_road}/${project_name}/pvalues.txt \
            --output-path=${output_road}/${project_name} \
            --output-name=${project_name}_plot.pdf 
            #--rows: File with a list of rows to plot, one per line [all available]
            #--columns: File with a list of columns to plot, one per line [all available]
            #--verbose / --quiet: Print or hide CellPhoneDB logs [verbose] 

$cellphonedb plot heatmap_plot \
             $input_road/${meta_file} \
             --pvalues-path=${output_road}/${project_name}/pvalues.txt \
             --output-path=${output_road}/${project_name} \
             --count-name=${project_name}_heatmap_count.pdf\
             --log-name=${project_name}_heatmap_log_count.pdf\
             --count-network-name=count_network.txt \
             --interaction-count-name=interactions_count.txt \
             --pvalue=0.05
#--verbose / --quiet: Print or hide cellphonedb logs [verbose]  
