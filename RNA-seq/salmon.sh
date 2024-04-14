
source /rsrch3/home/genomic_med/aksaw/anaconda3/bin/activate
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.8.0/salmon-1.8.0_linux_x86_64.tar.gz

tar -xf  salmon-1.8.0_linux_x86_64.tar.gz

# reference transcript latest release, 

mkdir gencode_refrence_transcript

cd gencode_refrence_transcript

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M28/gencode.vM28.transcripts.fa.gz

cd ..

mkdir mouse_reference_transcript_index

salmon-1.8.0_linux_x86_64/bin/salmon  index --gencode -t gencode_refrence_transcript/gencode.vM28.transcripts.fa.gz    -i mouse_reference_transcript_index

salmon-1.8.0_linux_x86_64/bin/salmon  index --gencode -t gencode_refrence_transcript/gencode.v39.transcripts.fa.gz    -i human_reference_transcript_index


cd  result

../salmon-1.8.0_linux_x86_64/bin/salmon  quant  -i ../human_reference_transcript_index  -l A  -1 ../fastq_HML_RNA_051022/21116XR-03-13_S40_L001_R1_001.fastq.gz  -2 ../fastq_HML_RNA_051022/21116XR-03-13_S40_L001_R2_001.fastq.gz  -p 10 --gcBias --numGibbsSamples 20  --validateMappings -o quants/S40_L001_quant


../salmon-1.8.0_linux_x86_64/bin/salmon  quant  -i ../human_reference_transcript_index  -l A -1 ../fastq_HML_RNA_051022/21116XR-03-14_S64_L001_R1_001.fastq.gz  -2 ../fastq_HML_RNA_051022/21116XR-03-14_S64_L001_R2_001.fastq.gz  -p 10 --gcBias --numGibbsSamples 20  --validateMappings -o quants/S64_L001_quant


../salmon-1.8.0_linux_x86_64/bin/salmon  quant  -i ../human_reference_transcript_index  -l A -1 ../fastq_HML_RNA_051022/21116XR-03-15_S0_L001_R1_001.fastq.gz  -2 ../fastq_HML_RNA_051022/21116XR-03-15_S0_L001_R2_001.fastq.gz  -p 10 --gcBias --numGibbsSamples 20  --validateMappings -o quants/S0_L001_quant

../salmon-1.8.0_linux_x86_64/bin/salmon  quant  -i ../human_reference_transcript_index  -l A -1 ../fastq_HML_RNA_051022/21116XR-03-16_S52_L001_R1_001.fastq.gz  -2 ../fastq_HML_RNA_051022/21116XR-03-16_S52_L001_R2_001.fastq.gz  -p 10 --gcBias --numGibbsSamples 20  --validateMappings -o quants/S52_L001_quant

../salmon-1.8.0_linux_x86_64/bin/salmon  quant  -i ../human_reference_transcript_index  -l A -1 ../fastq_HML_RNA_051022/21116XR-03-17_S65_L001_R1_001.fastq.gz  -2 ../fastq_HML_RNA_051022/21116XR-03-17_S65_L001_R2_001.fastq.gz  -p 10 --gcBias --numGibbsSamples 20  --validateMappings -o quants/S65_L001_quant


../salmon-1.8.0_linux_x86_64/bin/salmon  quant  -i ../human_reference_transcript_index  -l A -1 ../fastq_HML_RNA_051022/21116XR-03-18V1_S110_L001_R1_001.fastq.gz  -2 ../fastq_HML_RNA_051022/21116XR-03-18V1_S110_L001_R2_001.fastq.gz  -p 10 --gcBias --numGibbsSamples 20  --validateMappings -o quants/S110_L001_quant


Rscript differential_analysis.R
