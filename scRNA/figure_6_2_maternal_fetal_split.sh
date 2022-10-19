
#BAM_FILE=A1_target_region.sorted.bam
#filter_file=A1_fetal_barcode.txt
input_bam_road=/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure6/GEX_Abortion_SNP_bam
filter_file_road=/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure6/GEX_Abortion_SNP_bam/barcode_split
out_road=/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure6/GEX_Abortion_SNP_bam/fetal_maternal_split_bam

#list="A1 A3 A4 A8 N1 N2 N3 N4 N5"

list="N1"
for sample in $list;
do
##for fetal
# Save the header lines
samtools view -H $input_bam_road/${sample}_target_region.sorted.bam > $out_road/${sample}_SAM_header

##for fetal
# Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
samtools view  $input_bam_road/${sample}_target_region.sorted.bam | LC_ALL=C grep -F -f $filter_file_road/${sample}_fetal_barcode.txt >  $out_road/${sample}_fetal_filtered_SAM_body
# Combine header and body
cat  $out_road/${sample}_SAM_header  $out_road/${sample}_fetal_filtered_SAM_body >  $out_road/${sample}_fetal_filtered.sam
# Convert filtered.sam to BAM format
samtools view -b  $out_road/${sample}_fetal_filtered.sam >  $out_road/${sample}_fetal_filtered.bam
samtools sort  $out_road/${sample}_fetal_filtered.bam  -o   $out_road/${sample}_fetal_filtered_sort.bam
samtools index $out_road/${sample}_fetal_filtered_sort.bam

##for maternal
# Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
samtools view  $input_bam_road/${sample}_target_region.sorted.bam | LC_ALL=C grep -F -f $filter_file_road/${sample}_maternal_barcode.txt >  $out_road/${sample}_maternal_filtered_SAM_body
# Combine header and body
cat  $out_road/${sample}_SAM_header  $out_road/${sample}_maternal_filtered_SAM_body >  $out_road/${sample}_maternal_filtered.sam
# Convert filtered.sam to BAM format
samtools view -b  $out_road/${sample}_maternal_filtered.sam >  $out_road/${sample}_maternal_filtered.bam
samtools sort  $out_road/${sample}_maternal_filtered.bam  -o   $out_road/${sample}_maternal_filtered_sort.bam
samtools index $out_road/${sample}_maternal_filtered_sort.bam

rm $out_road/${sample}_*_filtered.bam $out_road/${sample}_SAM_header $out_road/${sample}_*_filtered_SAM_body $out_road/${sample}_*_filtered.sam
done
#ref source https://kb.10xgenomics.com/hc/en-us/articles/360022448251-Is-there-way-to-filter-the-BAM-file-produced-by-10x-pipelines-with-a-list-of-barcodes-
