zcat E250004601_L01_24052000063_1.fq.gz | sed 's/\// /' -  | gzip - > E250004601_L01_24052000063_new_1.fq.gz
zcat E250004601_L01_24052000063_2.fq.gz | sed 's/\// /' -  | gzip - > E250004601_L01_24052000063_new_2.fq.gz

mkdir raw
mv *_new_*fq.gz raw/
cd raw/

for i in ./*_1.fq.gz; do base=$(basename $i "_1.fq.gz");  umi_tools whitelist --stdin $i --bc-pattern='(?P<cell_1>.{8})(?P<discard_1>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=1}(?P<cell_2>.{6})(?P<discard_2>CCCATGATCGTCCGATCGTCGGCAGCGTCTCCACGC){e<=1}(?P<cell_3>.{6})(?P<umi_1>.{8}).*'  --extract-method=regex --set-cell-number 1000 --ed-above-threshold=correct  --error-correct-threshold=2 --log2stderr --knee-method=distance > ${base}_whitelist.txt --plot-prefix ${base} ; done
cd ..

mkdir 01_extract

cd raw/
for i in ./*_whitelist.txt; do base=$(basename $i "_whitelist.txt");  umi_tools extract --bc-pattern='(?P<cell_1>.{8})(?P<discard_1>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=1}(?P<cell_2>.{6})(?P<discard_2>CCCATGATCGTCCGATCGTCGGCAGCGTCTCCACGC){e<=1}(?P<cell_3>.{6})(?P<umi_1>.{8}).*' --stdin ${base}_1.fq.gz  --stdout ../01_extract/${base}_1.extract.fq.gz  --read2-in  ${base}_2.fq.gz  --read2-out ../01_extract/${base}_2.extract.fq.gz  --error-correct-cell  --extract-method=regex  --whitelist ${base}_whitelist.txt; done
cd ../01_extract/

for i in `ls *_2.extract.fq.gz`; do base=$(basename $i "_2.extract.fq.gz"); cutadapt -u 52 -o ${base}.trimmed.2.extract.fq.gz ${base}_2.extract.fq.gz; done

for i in `ls *_2.extract.fq.gz`; do base=$(basename $i "_2.extract.fq.gz"); cutadapt -q 20 -O 10 -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -B CTGTCTCTTATACACATCTGACGCTGCCGACGA -m 10 --max-n 0.1 --trim-n -o ${base}.cut.1.fq.gz -p ${base}.cut.2.fq.gz ${base}_1.extract.fq.gz ${base}.trimmed.2.extract.fq.gz; done

for i in `ls *_1.extract.fq.gz`; do base=$(basename $i "_1.extract.fq.gz"); zcat ${base}.cut.1.fq.gz |awk '{if(NR%4 == 1){print ">" substr($0, 2)}}{if(NR%4 == 2){print}}' - > ${base}_1.extract.fa; done

for i in `ls *_2.extract.fq.gz`; do base=$(basename $i "_2.extract.fq.gz"); zcat ${base}.cut.2.fq.gz |awk '{if(NR%4 == 1){print ">" substr($0, 2)}}{if(NR%4 == 2){print}}' - > ${base}_2.extract.fa; done
rm *trimmed*
cd ..
mkdir 02_align
for i in ./01_extract/*_2.extract.fa; do base=$(basename $i "_2.extract.fa"); bowtie2 -p 2 -f  -N 1 --very-sensitive-local --no-unal  -x /media/helab/data1/00_public/database/genomes/hg19_bowtie2/hg19 -1 ./01_extract/${base}_1.extract.fa -2 ./01_extract/${base}_2.extract.fa -S 02_align/${base}.hg19.sam 2> 02_align/${base}.hg19.align.log; done

cd ./01_extract
rm *fa
cd ..
mkdir 03_bam
for i in 02_align/*.sam; do base=$(basename $i ".sam"); samtools view -hbS -q 20 ./02_align/${base}.sam | samtools sort -T ${base} - > 03_bam/${base}.bam; done

cd 02_align
rm *sam
cd ..

mkdir 04_rmdup
for i in 03_bam/*.bam; do base=$(basename $i ".bam"); java -jar /media/helab/data1/00_public/software/picard-tools-2.2.4/picard.jar  MarkDuplicates REMOVE_DUPLICATES=true I=./03_bam/${base}.bam o=./04_rmdup/${base}_rmdup.bam M=./04_rmdup/${base}_rmdup_picard.txt; done

for i in ./04_rmdup/*_rmdup.bam; do base=$(basename $i "_rmdup.bam"); samtools view ./04_rmdup/${base}_rmdup.bam | sed 's/_/\tCB:Z:/'  - | sed 's/_/\tTB:Z:/' - |awk '{print$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$2"\t"$3}' -> ./04_rmdup/${base}_CB_TB.bam; done

for i in ./04_rmdup/*_CB_TB.bam; do base=$(basename $i "_CB_TB.bam"); samtools view ./04_rmdup/${base}_rmdup.bam  -H | cat - ./04_rmdup/${base}_CB_TB.bam | samtools sort -t CB - -o ./04_rmdup/${base}_CB_TB_sorted_tags.bam; done

cd 04_rmdup/
rm *_CB_TB.bam
cd ..

mkdir 05_split











