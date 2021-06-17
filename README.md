# Bulk-RNAseq-pipline(hisat2-StringTie-DEseq2)
## Download data via SRAToolkit
```
cd /export2/liuhw/Bulk_RNA_seq/workflow_test/
mkdir SRP200940
cd ./SRP200940
/export/bioinfo-team/home/liuhw/software/SRAToolkit/sratoolkit.2.10.9-ubuntu64/bin/prefetch \
    --option-file /export2/liuhw/Bulk_RNA_seq/workflow_test/SRR_Acc_List.txt
```

## Move srr file to a new file SRR_all
```
cd ./SRR_all
```

## SRA file to fastq 
```
mkdir fastq
##nohup bash rna.sh  > srr2fastq.log 2>&1 &

file_path=/export2/liuhw/Bulk_RNA_seq/workflow_test/SRP200940/
indir=${file_path}/SRR_all
outdir=${file_path}/fastq

for file in `ls $indir | grep .sra`;
do
 a=${file%.sra*};
 echo ----------
 echo ${a}
 /export/bioinfo-team/home/liuhw/software/SRAToolkit/sratoolkit.2.10.9-ubuntu64/bin/fastq-dump $indir/${a}.sra --gzip --split-files -O $outdir/
done
```
* #### 只有一个fastq为Single end 
* #### 两个fastq为Paired end 

## Trimmomatic remove low quality sequence and adapter
```
mkdir clean_fastq
##nohup bash rna.sh  > clean.log 2>&1 &

file_path=/export2/liuhw/Bulk_RNA_seq/workflow_test/SRP200940/
indir=${file_path}/fastq
outdir=${file_path}/clean_fastq
 
for file in `ls $indir | grep 1.fastq.gz`;
do
 a=${file%_1.fastq.gz*}; 
 echo ----------
 echo ${a};
 java -jar /export2/liuhw/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 $indir/${a}_1.fastq.gz $outdir/${a}_1.clean.fastq.gz  ILLUMINACLIP:/export2/liuhw/software/Trimmomatic-0.39/adapters/TruSeq2-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:20 TRAILING:20 MINLEN:75 & 
done
```
* **Single End example**: 
```
java -jar trimmomatic-0.35.jar PE -phred33 input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
* **Paired End example**: 
```
java -jar trimmomatic-0.35.jar PE -phred33 input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
* Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
* Remove leading low quality or N bases (below quality 3) (LEADING:3)
* Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
* Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
* Drop reads below the 36 bases long (MINLEN:36)

## fastqc
```
# conda activate macs2
##nohup bash rna.sh  > fastq.log 2>&1 &

fastqc=/export2/liuhw/software/FastQC/fastqc 
 
cd ./clean_fastq
$fastqc *fastq.gz
multiqc .
```
* multiqc整合一下结果

## hisat2 and samtool
```
mkdir bam
##nohup bash rna.sh  > hista2.log 2>&1 &

file_path=/export2/liuhw/Bulk_RNA_seq/workflow_test/SRP200940/
hisat=/Shared_Software/Genomic/hisat2-2.2.1/hisat2
GRCh38=/Shared_Software/ref_genome/hisat2_grch38_tran/grch38_tran/genome_tran

outdir=${file_path}/bam
indir=${file_path}/clean_fastq
 
for file in `ls $indir | grep 1.clean.fastq.gz`;
do
 a=${file%_1.clean.fastq.gz*};
 echo -------------;
 echo ${a};
 $hisat -p 60 --dta -x $GRCh38 -U $indir/${a}_1.clean.fastq.gz | samtools view -Sbh - > $outdir/${a}.bam
done
```
* **Single End example**: 
```
hisat2 -f -x genome_index -U reads_1.fastq -S eg1.sam
```
* **Paired End example**: 
```
hisat -f -x genome_index -1 reads_1.fastq -2 reads_2.fastq -S eg2.sam
```
* 用samtool把sam转为bam

## sort bam
```
mkdir sorted_bam
##nohup bash rna.sh  > sorted_bam.log 2>&1 &

file_path=/export2/liuhw/Bulk_RNA_seq/workflow_test/SRP200940/
 
indir=${file_path}/bam
outdir=${file_path}/sorted_bam
 
for file in `ls $indir | grep .bam`;
do
 a=${file%.bam*};
 echo ----------
 echo ${a}
 date
 samtools sort -@ 30 -o $outdir/${a}.sort.bam $indir/${a}.bam
done
```

## StringTie
```
mkdir stringtie_gtf
mkdir ballgown
##nohup bash rna.sh  > StringTie.log 2>&1 &

file_path=/export2/liuhw/Bulk_RNA_seq/workflow_test/SRP200940/
stringtie=/Shared_Software/Genomic/stringtie-2.1.4.Linux_x86_64/stringtie
GRCh38_gtf=/Shared_Software/ref_genome/refdata-gex-GRCh38-2020-A/genes/genes.gtf

outdir=${file_path}/stringtie_gtf
indir=${file_path}/sorted_bam
ballgown_dir=${file_path}/ballgown

for file in `ls $indir | grep .sort.bam`;
do
 a=${file%.sort.bam*};
 echo -------------;
 echo ${a};
 date
 $stringtie -p 50 -G $GRCh38_gtf -o $outdir/${a}.gtf $indir/${a}.sort.bam
done
 
echo ----------
echo "merge gtf files"
date
cd $outdir
$stringtie --merge -p 50 -G $GRCh38_gtf  -o $file_path/merged.gtf `ls $outdir | grep .gtf`
echo "merge done"

echo ----------
echo "counting transcript"
for file in `ls $indir | grep .sort.bam`;
do
 a=${file%.sort.bam*};
 echo -------------;
 echo ${a};
 date
 $stringtie -e -B -p 50 -G $file_path/merged.gtf -o $ballgown_dir/${a}/${a}.gtf $indir/${a}.sort.bam
done
```

## ballgown to DEseq2
```
##nohup bash rna.sh  > prepDE.log 2>&1 &
prepDE=/Shared_Software/Genomic/stringtie-2.1.4.Linux_x86_64/prepDE.py

python2 $prepDE -i gtf_path.txt -g gene_count_matrix.csv -t transcript_count_matrix.csv -l 150
```

