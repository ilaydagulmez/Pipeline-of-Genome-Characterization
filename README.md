# Our current pipeline of genome characterization of polyploid plants

# Step 1: Quality Check

Use FastQC to quality check. If filter needed pass to step 2.

# Step 2: Filtering and Correction

Use SOAPnuke (or fastp) to filtering and fastp to correction.

SOAPnuke filter -n 0.01 -l 35 -q 0.1 -1 1.fastq.gz -C clean_1.fastq.gz -2 2.fastq.gz -D clean_2.fastq.gz --polyX 65 --minReadLen 65

fastp --in1 clean_1.fastq --out1 FASTP_1.fastq.gz --in2 clean_2.fastq --out2 FASTP_2.fastq.gz --correction

# Step 3: De novo Assembly 

Use Karect before assembly step the second error correction

karect -correct -threads=16 -matchtype=hamming -celltype=diploid -inputfile=FASTP_1.fastq -inputfile=FASTP_2.fastq

Then use MEGAHIT and SOAPdenovo2 to de novo assembly. If you don't want to use optimum kmer (21), use Kmergenie for identify best kmer.

megahit --no-mercy -1 karect_FASTP_1.fastq -2 karect_FASTP_2.fastq

SOAPdenovo-fusion -D -s config.txt -p 64 -K 21 (or identify) -g k21_1 -c ./megahit_out/final.contigs.fa
SOAPdenovo-127mer map -s config.txt -p 64 -g k21_1
SOAPdenovo-127mer scaff -p 64 -g k21_1 -F

##output file = k21_1.scafSeq

After assembly step use Redundans to eliminate redundancy and potentially generate additional scaffolding:

redundans.py -v -i karect_FASTP_1.fastq karect_FASTP_2.fastq -f k21_1.scafSeq  --nogapclosing -o redundans

##output file = scaffolds.reduced.fa


# 3 times polishing

Use Racon to polish your assembly. This tool has 3 steps and 3 times repeat. 

First, cat your 2 reads: 

cat karect_FASTP_1.fastq  karect_FASTP_2.fastq > karect_FASTP_1_2.fastq 

Then start the pipeline:

# step1 
python racon_preprocess.py karect_FASTP_1_2.fastq > final.racon.fastq

# step2 
minimap2 -ax sr scaffolds.reduced.fa final.racon.fastq > racon_aln.sam

# step3 
racon -m 8 -x -6 -g -8 -w 500 -t 64 final.racon.fastq racon_aln.sam scaffolds.reduced.fa > polished1.fasta

After getting polished1.fasta, you need to use as input from # step2:

# second time polishing

minimap2 -ax sr polished1.fa final.racon.fastq > 2_racon_aln.sam
racon -m 8 -x -6 -g -8 -w 500 -t 64 final.racon.fastq 2_racon_aln.sam polished1 > polished2.fasta

Again, after getting polished2.fasta, you need to use as input from # step2:

# last time polishing

minimap2 -ax sr polished2.fa final.racon.fastq > 3_racon_aln.sam
racon -m 8 -x -6 -g -8 -w 500 -t 64 final.racon.fastq 3_racon_aln.sam polished2 > polished3.fasta

Finally use RagTag to correct and scaffold misassemblies against as reference genome:

ragtag.py scaffold reference.fasta polished3.fasta

##output file = ragtag.scaffold.fasta


# Step 4: Assessment of Assemblies

Use BUSCO to assessment of assembly and quality:

busco -i ragtag.scaffold.fasta -l viridiplantae_odb10 -m genome -o ragtag_busco

To plot BUSCO:

python3 generate_plot.py -wd ragtag_busco 


# Step 5: Identify Transposable Elements

Use RepeatMasker if you use default TE librarys but we created and recommend own genus library:

BuildDatabase -name ragtag.scaffold.fasta
RepeatModeler -database ragtag.scaffold.fasta -threads 36 -LTRStruct > out.log
RepeatMasker -e rmblast -x -pa 36 -gff -lib consensi.fa.classified ragtag.scaffold.fasta -dir ./final_directory


# Step 6: Annotation and Gene Prediction

Use AUGUSTUS to genome annotation:

augustus --strand=both --species=arabidopsis --protein=on --cds=on --gff3=on --codingseq=on ragtag.scaffold.fasta > ragtag.scaffold.fasta

To get protein sequence:

getAnnoFasta.pl ragtag_polished3.gff







