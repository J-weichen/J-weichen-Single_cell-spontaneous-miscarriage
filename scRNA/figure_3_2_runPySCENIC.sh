## inputs
road=$1 
#/mnt/data/chenwei/jianghai/3.seurat_result/3.SCENIC/loom_result
sn=$2
#ALLcell_DEGs_selected_0719

f_loom_path_scenic=$road/s1_${sn}.loom

## outputs
grn_output=$road/s2_${sn}.adj.tsv
ctx_output=$road/s2_${sn}.reg.tsv
f_pyscenic_output=$road/s2_${sn}.pyscenic.loom

## reference
f_tfs=/mnt/data/chenwei/10X_placenta/10x_data/191125-merge_nonormalization_add/SCENIC/human_hg38/hs_hgnc_tfs.txt
f_motif_path=/mnt/data/chenwei/10X_placenta/10x_data/191125-merge_nonormalization_add/SCENIC/human_hg38/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
#f_db_names=/home/chenwei/10x_data/191125-merge_nonormalization_add/SCENIC/human_hg38/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
#f_db_names=/home/chenwei/10x_data/191125-merge_nonormalization_add/SCENIC/human_hg38/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather

f_db_names=`find /mnt/data/chenwei/10X_placenta/10x_data/191125-merge_nonormalization_add/SCENIC/human_hg38/ -name "hg38*.feather"`
###step1  Co-expression network 

arboreto_with_multiprocessing.py $f_loom_path_scenic $f_tfs --method grnboost2 --output $grn_output --num_workers 50 --seed 19921010

pyscenic ctx $grn_output $f_db_names --annotations_fname $f_motif_path --expression_mtx_fname $f_loom_path_scenic --output $ctx_output --num_workers 20
## Note:
## The reason why I didn't use `--mask_dropouts` parameters:
## In R version of SCENIC, mask dropouts is False (default when pySCENIC>0.9.6, current version: 0.10.1). 
## Here we match the behavior of R version.
## 10G for one threed
# ## set TMPDIR to current path, in case of no enough disk space on /tmp/
export TMPDIR=/mnt/data/chenwei/gongchen/4.SCENIC/temp/

pyscenic aucell  $f_loom_path_scenic $ctx_output --output $f_pyscenic_output --num_workers 10 --seed 19921010

cp  $road/s2_${sn}.reg.txt  $road/s2_${sn}.reg.tsv

#threads=30 ; min_regulon_size=5
sn2=$road/s3_${sn}
python /mnt/data/chenwei/gongchen/0.script/figure_three_step_3_postSCENIC.py $f_pyscenic_output $ctx_output $sn2 5 5

cp $road/s2_${sn}.reg.tsv $road/s2_${sn}.reg.txt

mkdir $road/$sn
mv $road/*_${sn}* $road/$sn
