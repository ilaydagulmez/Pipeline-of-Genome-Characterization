# Our current pipeline of genome characterization of polyploid plants

# Step 1: Quality Check

Use FastQC (Andrew 2010) to quality check. If filter needed pass to step 2.


# Step 2: Filtering and Correction

Use SOAPnuke (Chen et al. 2018) to filtering and fastp (Chen et al, 2018) to correction.

SOAPnuke filter -n 0.01 -l 35 -q 0.1 -1 1.fastq.gz -C clean_1.fastq.gz -2 2.fastq.gz -D clean_2.fastq.gz --polyX 65 --minReadLen 65

fastp --in1 clean_1.fastq --out1 FASTP_1.fastq.gz --in2 clean_2.fastq --out2 FASTP_2.fastq.gz --correction

# Step 3: Finding Optimum k-mer, Genome Size and Genome Heterozygosity

Use Kmergenie (Chikhi and Medvedev, 2013) for identify optimum k-mer value. After getting optimum k-mer, use this for Jellyfish (Marcais and Kingsford, 2011) as input.

Use SGA (Simpson and Durbin, 2012) to estimate genome size.

sga preprocess --pe-mode 1 FASTP_1.fastq.gz FASTP_2.fastq.gz > genome.fastq

sga index -a ropebwt -t 8 --no-reverse genome.fastq

sga preqc -t 8 genome.fastq > P.vulgaris.preqc

sga-preqc-report.py P.vulgaris.preqc P.vulgaris.preqc

After that, merge two files:

cat FASTP_1.fastq.gz FASTP_2.fastq.gz > FASTP_1_2.fastq.gz

gunzip FASTP_1_2.fastq.gz

jellyfish count -C -m optimum k-mer -s 1000000000 -t 10 FASTP_1_2.fastq -o reads.jf
jellyfish histo -t 10 reads.jf > reads.histo

Now, reads.histo is also input for GenomeScope2.0 (http://genomescope.org/genomescope2.0/) Load this this file and get the genome heterozygosity results.


# Step 4: De novo Assembly 

Use Karect (Allam et al, 2015) before assembly step the second error correction

karect -correct -threads=16 -matchtype=hamming -celltype=diploid -inputfile=FASTP_1.fastq -inputfile=FASTP_2.fastq

Then use MEGAHIT and SOAPdenovo2 (Li et al, 2015) to de novo assembly. If you don't want to use optimum kmer (21), use Kmergenie for identify best kmer.

megahit --no-mercy -1 karect_FASTP_1.fastq -2 karect_FASTP_2.fastq

SOAPdenovo-fusion -D -s config.txt -p 64 -K 21 (or identify) -g k21_1 -c ./megahit_out/final.contigs.fa
SOAPdenovo-127mer map -s config.txt -p 64 -g k21_1
SOAPdenovo-127mer scaff -p 64 -g k21_1 -F

##output file = k21_1.scafSeq

After assembly step use Redundans (Pryszcz and Gabaldón, 2016) to eliminate redundancy and potentially generate additional scaffolding:

redundans.py -v -i karect_FASTP_1.fastq karect_FASTP_2.fastq -f k21_1.scafSeq  --nogapclosing -o redundans

##output file = scaffolds.reduced.fa


# 3 times polishing

Use Racon (Vaser et al. 2017) to polish your assembly. This tool has 3 steps and 3 times repeat. 

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

Finally use RagTag (Alonge et al, 2022) to correct and scaffold misassemblies against as reference genome:

ragtag.py scaffold reference.fasta polished3.fasta

##output file = ragtag.scaffold.fasta


# Step 5: Assessment of Assemblies

Use BUSCO (Simão et al. 2015) to assessment of assembly and quality:

busco -i ragtag.scaffold.fasta -l viridiplantae_odb10 -m genome -o ragtag_busco

To plot BUSCO:

python3 generate_plot.py -wd ragtag_busco 


# Step 6: Identify Transposable Elements

Use RepeatMasker (https://github.com/rmhubley/RepeatMasker) if you use default TE librarys but we created and recommend own genus library (Price et al. 2005)

BuildDatabase -name ragtag.scaffold.fasta
RepeatModeler -database ragtag.scaffold.fasta -threads 36 -LTRStruct > out.log
RepeatMasker -e rmblast -x -pa 36 -gff -lib consensi.fa.classified ragtag.scaffold.fasta -dir ./final_directory


# Step 7: Annotation and Gene Prediction

Use AUGUSTUS (Stanke et al. 2004) to genome annotation:

augustus --strand=both --species=arabidopsis --protein=on --cds=on --gff3=on --codingseq=on ragtag.scaffold.fasta > ragtag.scaffold.fasta

To get protein and CDS sequence:

getAnnoFasta.pl ragtag_polished3.gff
##output file = ragtag_polished3.aa and ragtag_polished3.codingseq


# Step 8: WGD and Distribution of Substitutions Per Synonymous Site (Ks) Analysis

Use wgd v2 (Chen et al, 2024) and CDS file which obtained from AUGUSTUS to identify Ks value. If you want to change your steps you can follow the wgd v2 github file.

wgd dmd ragtag_polished3.codingseq -o dmd

wgd ksd ragtag_polished3.codingseq.tsv ragtag_polished3.codingseq -o ksd

wgd syn -f transcript -a ID ragtag_polished3.codingseq.tsv ragtag_polished3.gff -ks ragtag_polished3.codingseq.tsv.ks.tsv --pathiadhore ./i-adhore -o syn

wgd peak --heuristic ragtag_polished3.codingseq.tsv.ks.tsv -ap ./syn/iadhore-out/anchorpoints.txt -sm ./syn/iadhore-out/segments.txt -le ./syn/iadhore-out/list_elements.txt -mp ./syn/iadhore-out/multiplicon_pairs.txt -n 1 4 -kc 3 -o wgd_peak

(Add outgroups CDS file for this step)

wgd dmd -f ragtag_polished3.codingseq -ap wgd_peak/AnchorKs_FindPeak/ragtag_polished3.codingseq.tsv.ks.tsv_95%CI_AP_for_dating_weighted_format.tsv -o wgd_dmd_ortho *add_all_outgroups_cds_file_here_respectively*

wgd focus --protcocdating --aamodel lg wgd_dmd_ortho/merge_focus_ap.tsv -sp dating_tree.nw -o wgd_dating -d mcmctree -ds 'burnin = 2000' -ds 'sampfreq = 1000' -ds 'nsample = 20000'  *add_all_outgroups_cds_file_here_respectively* ragtag_polished3.codingseq

NOTE: You should prepare dating_tree.nw file manually. For details, you can follow wgd v2 github file again.


# Step 9: Analysis of Duplicated Gene Pairs

Use doubletrouble R package (Almeida-Silva and Van de Peer 2024) and protein sequences which obtained from AUGUSTUS to classify gene duplication modes in one of five categories: segmental (SD), small-scale (SSD), dispersed (DD), proximal (PD), tandem (TD) duplications.

For the R package tool you can follow the steps here: https://www.bioconductor.org/packages/release/bioc/vignettes/doubletrouble/inst/doc/doubletrouble_vignette.html


# Step 10: Chloroplast Genome Characterization (optional)

Use GetOrganelle (Jin et al, 2020) to assembly of raw reads:

get_organelle_from_reads.py -1 karect_FASTP_1.fastq  -2 karect_FASTP_2.fastq  -o plastome_output -R 15 -k 21,45,65,85,105 -F embplant_pt

Use GeSeq software (https://chlorobox.mpimp-golm.mpg.de/geseq.html) to annotation. Then use OrganellarGenomeDRAW (OGDRAW) (https://chlorobox.mpimp-golm.mpg.de/OGDraw.html)  to draw circular maps.

Also use CPGAVAS2 (http://47.96.249.172:16019/analyzer/home) to identify codon usage. Finally use MISA (http://pgrc.ipk-gatersleben.de/misa/) to detect SSRs.

 














