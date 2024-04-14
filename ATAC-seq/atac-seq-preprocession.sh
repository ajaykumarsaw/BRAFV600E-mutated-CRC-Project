
for f1 in *R1.fastq.gz
do

fastqc  $f1  ${f1%%R1.fastq.gz}R2.fastq.gz

trim_galore --paired --fastqc -o Trim_glore/   $f1  ${f1%%R1.fastq.gz}R2.fastq.gz

bowtie2 --very-sensitive -X 2000 -x /rsrch3/home/genomic_med/aksaw/ajay/UCSC_REF_GENOME_HG19_bowtie2/bowtie2/hg19_index  -1 Trim_glore/${f1%%R1.fastq.gz}R1_val_1.fq.gz  -2  Trim_glore/${f1%%R1.fastq.gz}R2_val_2.fq.gz  -p 10  | samtools view -u - | samtools sort -  >  ${f1%%_*}.bam


# Mitochondrial reads removal 


samtools view -h ${f1%%_*}.bam  |  python3 atacseq/removeChrom.py - - chrM  |  samtools view -b -  > ${f1%%_*}.rmChrM.bam


samtools sort ${f1%%_*}.rmChrM.bam   -o  ${f1%%_*}.rmChrM.sorted.bam

samtools index ${f1%%_*}.rmChrM.sorted.bam

# remove PCR duplicates

java -jar picard.jar MarkDuplicates I= ${f1%%_*}.rmChrM.sorted.bam  O= rem_dup_${f1%%_*}.rmChrM.sorted.bam  M=${f1%%_*}_dups.txt REMOVE_DUPLICATES=true

samtools view -b -q 10 rem_dup_${f1%%_*}.rmChrM.sorted.bam >  rmMulti_rem_dup_${f1%%_*}.rmChrM.sorted.bam

samtools view -h rmMulti_rem_dup_${f1%%_*}.rmChrM.sorted.bam  |  python3 atacseq/SAMtoBED.py  -i -  -o ${f1%%_*}.bed  -x  -v

# Peak calling

macs2 callpeak  -t ${f1%%_*}.bed  -f BEDPE  -n ${f1%%_*}  -g hs  --keep-dup all

done

