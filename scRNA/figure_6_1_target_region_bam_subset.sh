bam_road=/mnt/data/chenwei/gongchen/0.script
out_road=/mnt/data/chenwei/gongchen/manuscript/manu_table/major_figure6/GEX_Abortion_SNP_bam
list="A1 A3 A4 A8 N1 N2 N3 N4 N5"
for sample in $list;
do
samtools view -hb -L  /mnt/data/chenwei/gongchen/manuscript/manu_table/Abortion_SNP_308_5kb_hg38_sort_new.bed  $bam_road/${sample}_out/outs/${sample}_possorted_genome.bam > $out_road/${sample}_target_region.bam
cd  $out_road
samtools sort ${sample}_target_region.bam  -o  ${sample}_target_region.sorted.bam
samtools index ${sample}_target_region.sorted.bam
done

###COMMAND TO VIEW THE TAG IN 10X BAM ： https://cloud.tencent.com/developer/article/1677573
#samtools view outs/possorted_bam.bam | head | perl -nle '@reads=split(/\t/,$_); if (m/CB:Z:([^\t\n]+)\t/) { print "$reads[2]\t","$reads[3]\t","$reads[3]\t","$reads[0]_",$1; }
##or 
#samtools view outs/possorted_bam.bam | head | awk -v OFS="\t" -F" " '{ print $3,$4,$4,$1";"$19}'


