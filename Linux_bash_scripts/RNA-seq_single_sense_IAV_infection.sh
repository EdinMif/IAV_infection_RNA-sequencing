##########################################################################
# RNA-seq Time Course: IAV infected mice time course experiment  		 #
# single-end reads.                                                      #
#       --- Linux bioinformatics workflow for known sense genes ---      #
##########################################################################
# Based on the pipeline created by Nalpas, N.C. (2014) and Correia, C.N
# DOI badge: http://dx.doi.org/10.5281/zenodo.12474
# Author of current version (4.0.1): Mifsud, E.J.
# DOI badge of current version:
# Last updated on: 08/03/2017

########################################
# File Transfer and change permissions #
########################################

# create directory for raw files
mkdir -p $HOME/IAV/rawfiles
cd !$

# All files were downloaded from the CZC server in 2017. 
# Files were transfered to the supercomputer in a folder reflecting time point
# post infection
   
# Change directory permissions to read and execute only:
for file in `find $HOME/IAV/rawfiles/D3 -name *.fastq`; \
do chmod 555 $file; done

# compress fastq files
for file in `find $HOME/IAV/rawfiles/D3 -name *.fastq`; \
do gzip -9 $file; done


# File names correspond to:
# labCode_TreatmentGroup_Timepoint_fastq.gz


###########################################
# FastQC quality check of raw FASTQ files #
###########################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir -p $HOME/IAV/quality_check/pre-filtering/chip2 
cd !$

# To create a job on the CZC supercomputer we must use the code bsub -q C
# Run FastQC in one file to see if it's working well:


bsub -q C fastqc -o /home/ejmifsud/IAV/quality_check/pre-filtering/ \
--noextract --nogroup -t 1 \
/home/ejmifsud/IAV/rawfiles/D3/R_2017_02_02_16_06_09_sn247770192_sn247770192-4-170202_EJM_D3_pi__/\
IonXpressRNA_001.R_2017_02_02_16_06_09_sn247770192_sn247770192-4-170202_EJM_D3_pi__.fastq.gz 

### Moved this folder to my laptop using Filezilla.
### and checked the HTML report. It worked fine.

# Create a bash script to perform FastQC quality check on all fastq.gz files:
cd $HOME/IAV/quality_check/pre-filtering

for file in `find $HOME/IAV/rawfiles/D3/ \
-name *fastq.gz`; do echo "fastqc --noextract --nogroup -t 5 \
-o $HOME/IAV/quality_check/pre-filtering $file" \
>> fastqc.sh; done
# count the number of lins in the file to make sure the correct files have been created
wc -l fastqc.sh 
nano fastqc.sh
#!/bin/sh
#BSUB -J fastqc_job
#BSUB -n 5
bsub -q C < fastqc.sh

# Deleted all the HTML files:
rm -r *.html

# because we process the samples using the queue to process the job the standerr has 
#all the information that would be shown as in the nohub folder
# This information contains the fastQC run information and must be looked at to ensure all
# files have successfully been run
# Check if all the files were processed:

more 442391.standout | grep 'Analysis complete' | wc -l
  
# Files were transfered to computer using Filezilla 

# Check all output from FastQC:
mkdir $HOME/IAV/quality_check/pre-filtering/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/IAV/quality_check/pre-filtering/tmp; \
done;

for file in \
`find $HOME/IAV/quality_check/pre-filtering/tmp \
-name summary.txt`; do more $file >> reports_pre-filtering.txt; done

for file in \
`find $HOME/IAV/quality_check/pre-filtering/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_pre-filtering.txt; \
done

# Remove temporary folder and its files:
rm -r $HOME/IAV/quality_check/pre-filtering/tmp

### After checking the number of sequences it was apparent there was an issue with 
### Chip-2 and the file transfer was corrupted. Therefore the process of fastQC analysis 
### for Chip-2 had to be repeated
# Files were downloaded from ion torrent server and transfered from my macbook air 
# to the CZC supercomputer using fileZilla

# Create a new folder for the new raw files 
mkdir  -p $HOME/IAV/rawfiles/D3chip2
cd !$

# Change directory permissions to read and execute only:
for file in `find $HOME/IAV/rawfiles/D3chip2 -name *.fastq`; \
do chmod 555 $file; done

# compress fastq files

for file in `find $HOME/IAV/rawfiles/D3chip2/R_2017_02_14_EJM_D3_pi_take_2_ -name *.fastq`; \
do echo 'gzip -9 -o $HOME/IAV/rawfiles/D3chip2 $file'; >> compress.sh; done

wc -l compress.sh 
nano compress.sh
#!/bin/sh
#BSUB -J compress_job
#BSUB -n 3
bsub -q C < compress.sh 
# all of the above failed so files were zipped overnight using the following loop

for file in `find $HOME/IAV/rawfiles/D3chip2 -name *.fastq`; \
do gzip -9 $file; done


########################
# Fastqc quality check #
########################

cd $HOME/IAV/quality_check/pre-filtering

for file in `find $HOME/IAV/rawfiles/D3chip2/ \
-name *fastq.gz`; do echo "fastqc --noextract --nogroup -t 3 \
-o $HOME/IAV/quality_check/pre-filtering $file" \
>> fastqcD3c2.sh; done

wc -l fastqc.sh 
nano fastqc.sh
#!/bin/sh
#BSUB -J fastqc_job
#BSUB -n 3
bsub -q C < fastqcD3c2.sh

#check the standout

more 443622.stdout | grep 'Analysis complete' | wc -l

# remove html files 
rm -r *.html

# Check all output from FastQC:
mkdir $HOME/IAV/quality_check/pre-filtering/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/IAV/quality_check/pre-filtering/tmp; \
done;

for file in \
`find $HOME/IAV/quality_check/pre-filtering/tmp \
-name summary.txt`; do more $file >> reports_pre-filtering.txt; done

for file in \
`find $HOME/IAV/quality_check/pre-filtering/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_pre-filtering.txt; \
done

# Remove temporary folder and its files:
rm -r $HOME/IAV/quality_check/pre-filtering/tmp

# Once files were checked the contents of D3chip D3chip2/R_2017_02_14_EJM_D3_pi_take_2
# fastqc files were moved to D3 folder using Filezilla and the corrupt files were deleted 


##################################################################
# Adapter-contamination and quality filtering of raw FASTQ files #
##################################################################

# Required software is cutadapts (version 1.12). More information can be found
# here: http://cutadapt.readthedocs.io/en/stable/guide.html

# Create a working directory for filtered reads:
mkdir $HOME/IAV/fastq_sequence/
cd $HOME/IAV/fastq_sequence/

# Run cutadapt in one file to see if it's working well and determine the parameters:

bsub -q C cutadapt -a ATCACCGACTGCCCATAGAGAGGCTGAGAC -g CTAAGGTAAC -g CCTCTCTATGGGCAGTCGGTGAT \
-e 0.3 -m 50 --quality-base 20 --discard-trimmed -o /home/ejmifsud/IAV/fastq_sequence.fastq \
/home/ejmifsud/IAV/rawfiles/D3/R_2017_02_02_16_06_09_sn247770192_sn247770192-4-170202_EJM_D3_pi__/\
IonXpressRNA_001.R_2017_02_02_16_06_09_sn247770192_sn247770192-4-170202_EJM_D3_pi__.fastq.gz

#check the fastq file 
bsub -q C fastqc -o /home/ejmifsud/IAV/fastq_sequence/ \
--noextract --nogroup -t 1 \
/home/ejmifsud/IAV/fastq_sequence.fastq

# Trial number 2 cutadapt
bsub -q C cutadapt -a ATCACCGACTGCCCATAGAGAGGCTGAGAC -g CCTCTCTATGGGCAGTCGGTGAT \
-m 50 --quality-base 20 --discard-trimmed \
-o /home/ejmifsud/IAV/fastq_sequence/1_PBS1_D3.fastq \
/home/ejmifsud/IAV/rawfiles/D3/R_2017_02_02_16_06_09_sn247770192_sn247770192-4-170202_EJM_D3_pi__/\
IonXpressRNA_001.R_2017_02_02_16_06_09_sn247770192_sn247770192-4-170202_EJM_D3_pi__.fastq.gz

bsub -q C fastqc -o /home/ejmifsud/IAV/fastq_sequence/ \
--noextract --nogroup -t 1 \
/home/ejmifsud/IAV/fastq_sequence/1_PBS1_D3.fastq

# Trial Number 3 cutadapt 
bsub -q C cutadapt -a ATCACCGACTGCCCATAGAGAGGCTGAGAC -g CCTCTCTATGGGCAGTCGGTGAT \
-m 50 -q 20,20 --discard-trimmed \
-o /home/ejmifsud/IAV/fastq_sequence/1_PBS1_D3.fastq \
/home/ejmifsud/IAV/rawfiles/D3/R_2017_02_02_16_06_09_sn247770192_sn247770192-4-170202_EJM_D3_pi__/\
IonXpressRNA_001.R_2017_02_02_16_06_09_sn247770192_sn247770192-4-170202_EJM_D3_pi__.fastq.gz

bsub -q C fastqc -o /home/ejmifsud/IAV/fastq_sequence/ \
--noextract --nogroup -t 1 \
/home/ejmifsud/IAV/fastq_sequence/1_PBS1_D3.fastq

# Trial Number 4 

bsub -q C cutadapt -a ATCACCGACTGCCCATAGAGAGGCTGAGAC -g CCTCTCTATGGGCAGTCGGTGAT \
-m 50 -q 24,24 --discard-trimmed \
-o /home/ejmifsud/IAV/fastq_sequence/1_PBS1_D3_2.fastq \
/home/ejmifsud/IAV/rawfiles/D3/R_2017_02_02_16_06_09_sn247770192_sn247770192-4-170202_EJM_D3_pi__/\
IonXpressRNA_001.R_2017_02_02_16_06_09_sn247770192_sn247770192-4-170202_EJM_D3_pi__.fastq.gz

bsub -q C fastqc -o /home/ejmifsud/IAV/fastq_sequence/ \
--noextract --nogroup -t 1 \
/home/ejmifsud/IAV/fastq_sequence/1_PBS1_D3_2.fastq


# Trial number 6 

bsub -q C cutadapt -q 24,24 -a ATCACCGACTGCCCATAGAGAGGCTGAGAC -g CCTCTCTATGGGCAGTCGGTGAT \
-n 5 -m 50 -M 240 --discard-trimmed \
-o /home/ejmifsud/IAV/fastq_sequence/1_PBS1_D3_4.fastq \
/home/ejmifsud/IAV/rawfiles/D3/R_2017_02_02_16_06_09_sn247770192_sn247770192-4-170202_EJM_D3_pi__/\
IonXpressRNA_001.R_2017_02_02_16_06_09_sn247770192_sn247770192-4-170202_EJM_D3_pi__.fastq.gz

bsub -q C fastqc -o /home/ejmifsud/IAV/fastq_sequence/ \
--noextract --nogroup -t 1 \
/home/ejmifsud/IAV/fastq_sequence/1_PBS1_D3_4.fastq



# Run cutadpt on all samples 

for file in `find $HOME/IAV/rawfiles/D3/ \
-name *fastq.gz`; \
do sample=`basename $file`; \
echo "cutadapt -q 24,24 -a ATCACCGACTGCCCATAGAGAGGCTGAGAC \
-g CCTCTCTATGGGCAGTCGGTGAT -n 5 -m 50 -M 240 --discard-trimmed \
-o $HOME/IAV/fastq_sequence/trimmed_${sample} $file" >> cutadapt.sh; done

wc -l cutadapt.sh 

nano cutadapt.sh
#!/bin/sh
#BSUB -J cutadapt_job
#BSUB -n 5

bsub -q C < cutadapt.sh


################################################
# FastQC quality check of filtered FASTQ files #
################################################
# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# Check all output from FastQC:
# make an output folder for Fastqc files 

mkdir $HOME/IAV/fastq_sequence/trimmed_fastq

for file in `find $HOME/IAV/fastq_sequence/ \
-name *fastq.gz`; do echo "fastqc --noextract --nogroup -t 5 \
-o $HOME/IAV/fastq_sequence/trimmed_fastq $file" \
>> trimmed_fastqc.sh; done

wc -l trimmed_fastqc.sh


nano trimmed_fastqc.sh
#!/bin/sh
#BSUB -J trimmed_fastq_job
#BSUB -n 5

bsub -q C < trimmed_fastqc.sh

more 446473.stdout | grep 'Analysis complete' |wc -l

# remove html files 
rm -r *.html


for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/IAV/fastq_sequence/trimmed_fastq/tmp; \
done;

for file in \
`find $HOME/IAV/fastq_sequence/trimmed_fastq/tmp \
-name summary.txt`; do more $file >> reports_pre-filtering.txt; done

for file in \
`find $HOME/IAV/fastq_sequence/trimmed_fastq/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_pre-filtering.txt; \
done

################################################################################
# Alignment of FASTQ files against the Mus_musculus reference genome with STAR #
################################################################################

# Required software is STAR 2.5.1b, consult manual/tutorial for details:
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

# Download Mus_musculus genome from Ensemble website 

mkdir -p $HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/source_file
#exit the internal node to enter the management node and download the file 
exit
cd $HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/source_file
wget ftp://ftp.ensembl.org/pub/release-87/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz

#log back into the internal node 

cd $HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/source_file
gunzip Mus_musculus.GRCm38.dna.toplevel.fa.gz


# Download Mus_musculus genome from Ensemble website 

mkdir -p $HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/source_file
#exit the internal node to enter the management node and download the file 
exit
cd $HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/source_file
wget ftp://ftp.ensembl.org/pub/release-87/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz

#log back into the internal node 

cd $HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/source_file
gunzip Mus_musculus.GRCm38.dna.toplevel.fa.gz

# For future reference this is the website to download the genome 
# ftp://ftp.ensembl.org/pub/release-87/gtf/mus_musculus/
# Download annotation file for Mus_musculus genome version GRC38 Annotation Release-87:
mkdir $HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/annotation_file

# exit internal node 
exit 
cd $HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/annotation_file
wget ftp://ftp.ensembl.org/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf.gz
 
# Log back into the Internal node 
 cd cd $HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/annotation_file
 gunzip Mus_musculus.GRCm38.87.gtf.gz
 
# For future reference this is the website to download annotation:
# ftp://ftp.ensembl.org/pub/release-87/fasta/mus_musculus/dna/ 
# Generate genome indexes files using annotations:
mkdir $HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/STAR-2.5.1b_index
cd $HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/STAR-2.5.1b_index

# NB: sjdboverhang should not exceed 100 as this does not help alignment and will 
#cause erros to occur when indexing the genome
bsub -q C STAR --runThreadN 20 --runMode genomeGenerate \
--genomeDir /$HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/STAR-2.5.1b_index \
-- genomeFastaFiles \
$HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/source_file/Mus_musculus.GRCm38.dna.toplevel.fa \
--sjdbGTFfile /$HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/annotation_file/Mus_musculus.GRCm38.87.gtf \
--sjdbOverhang 100 \
--outFileNamePrefix  \
$HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/STAR-2.5.1b_index


# Create and enter alignment working directory:
mkdir $HOME/IAV/STAR-2.5.1b_alignment
cd $HOME/IAV/STAR-2.5.1b_alignment

# Mapping reads from one FASTQ file to the indexed genome,
# to check if it works well:

bsub -q C STAR --runMode alignReads --runThreadN 1 --genomeLoad LoadAndRemove \
--genomeDir $HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/STAR-2.5.1b_index \
--readFilesIn \
$HOME/IAV/fastq_sequence/trimmed_IonXpressRNA_001.R_2017_02_02_16_06_09_sn247770192_sn247770192-4-170202_EJM_D3_pi__.fastq.gz \
--readFilesCommand gunzip -c --outFilterMultimapNmax 10 \
--outFilterMismatchNmax 10 --outFileNamePrefix ./trimmed_IonXpressRNA_ \
--outSAMtype BAM Unsorted --outReadsUnmapped Fastx 

# FASTQC of Trial 

bsub -q C fastqc -o /home/ejmifsud/IAV/STAR-2.5.1b_alignment/ \
--noextract --nogroup -t 1 \
/home/ejmifsud/IAV/fastq_sequence/1_PBS1_D3_4.fastq

# Script for alignment 

for file in `find $HOME/IAV/fastq_sequence/ \
-name *fastq.gz`; \
do sample=`basename $file | perl -p -e 's/_2017_.+\.fastq\.gz//'`; \
foldername=`basename $sample | perl -p -e 's/trimmed\_//'`; \
echo "mkdir $HOME/IAV/STAR-2.5.1b_alignment/$foldername; \
cd $HOME/IAV/STAR-2.5.1b_alignment/$foldername; \
STAR --runMode alignReads --runThreadN 20 --genomeLoad LoadAndRemove \
--genomeDir $HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/STAR-2.5.1b_index \
--readFilesIn $file --readFilesCommand gunzip -c \
--outFilterMultimapNmax 10 --outFilterMismatchNmax 10 \
--outFileNamePrefix ./$foldername --outSAMtype BAM Unsorted \
--outSAMattrIHstart 0 --outSAMattributes Standard --outReadsUnmapped Fastx" \
>> alignment.sh; \
done

nano alignment.sh
#!/bin/sh
#BSUB -J alignment_job
#BSUB -n 20

bsub -q C <  alignment.sh

# check the JOB#.stdout in the output file to ensure that the job was complete 
#Perl script Star report opener was forked from carol Coreia github account and then 
#transfered from laptop to the CZC supercomputer using Filezilla and placed 
# in a scripts folder 
# Merge all STAR log.final.out files into a single file:
for file in `find $HOME/IAV/STAR-2.5.1b_alignment \
-name *Log.final.out`; \
do perl /$HOME/Scripts/star_report_opener.pl -report $file; done;

# star_report_opener.pl report was transfered to the laptop using filezilla 


#############################################
# FastQC quality check of aligned BAM files #
#############################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and go to working directory:
mkdir $HOME/IAV/quality_check/post_alignment
cd $HOME/IAV/quality_check/post_alignment

# Create a bash script to perform FastQC quality check on aligned SAM files:
for file in `find $HOME/IAV/STAR-2.5.1b_alignment \
-name *.bam`; do echo "fastqc --noextract --nogroup -t 2 \
-o $HOME/IAV/quality_check/post_alignment $file" >> \
fastqc_aligned.sh; done;

ls 
wc -l fastqc_aligned.sh 

nano fastqc_aligned.sh
#!/bin/sh
#BSUB -J fastqc_aligned_job
#BSUB -n 10

bsub -q C < fastqc_aligned.sh
# Check the Job#.stdout to ensure that the job is complete
more 448741.stdout

# Delete all the HTML files:
rm -r *.html

# Check all output from FastQC:
mkdir $HOME/IAV/quality_check/post_alignment/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/IAV/quality_check/post_alignment/tmp; \
done

for file in \
`find $HOME/IAV/quality_check/post_alignment/tmp \
-name summary.txt`; do more $file >> reports_post-alignment.txt; done

for file in \
`find $HOME/IAV/quality_check/post_alignment/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_post_alignment.txt; \
done

# Check if all files were processed:
grep -c '##FastQC' basic_stats_post_alignment.txt
grep -c 'Basic Statistics' reports_post-alignment.txt

# Remove temporary folder:
rm -r tmp/

###################################################################
# Summarisation of gene counts with featureCounts for sense genes #
###################################################################

# Required package is featureCounts, which is part of Subread 1.5.0-p1 software,
# consult manual for details:
# http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf

# Create working directories:
mkdir $HOME/IAV/count_summerisation
cd  $HOME/IAV/count_summerisation

# Run featureCounts with one sample to check if it is working fine:
bsub -q C featureCounts -a \
$HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/annotation_file/Mus_musculus.GRCm38.87.gtf \
-R -s 1 -T 10 -t gene -g gene_id -o ./counts.text \
$HOME/IAV/STAR-2.5.1b_alignment/IonXpressRNA_001.R/IonXpressRNA_001.RAligned.out.bam

bsub -q C featureCounts -a \
$HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/annotation_file/Mus_musculus.GRCm38.87.gtf \
-R - s 1 -T 10 -t gene -g gene_id -o ./counts08.text \
$HOME/IAV/STAR-2.5.1b_alignment/IonXpressRNA_008.R/IonXpressRNA_008.RAligned.out.bam



# Create a bash script to run featureCounts on BAM file containing multihits and
# uniquely mapped reads using the reversely stranded parameter:
for file in `find $HOME/IAV/STAR-2.5.1b_alignment/ \
-name *Aligned.out.bam`; \
do sample=`basename $file | perl -p -e 's/Aligned\.out\.bam//'`; \
foldername=`basename $sample | perl -p -e 's/\.R\_//'`; \
echo "mkdir $HOME/IAV/count_summerisation/$foldername; \
cd $HOME/IAV/count_summerisation/$foldername; \
featureCounts -a \
$HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/annotation_file/Mus_musculus.GRCm38.87.gtf \
-R -s 1 -T 10 -t gene -g gene_id \
-o ${sample}_sense-counts.txt $file" >> sense_count.sh; done

nano sense_count.sh
#!/bin/sh
#BSUB -J feature_counts_job
#BSUB -n 10

bsub -q C < sense_count.sh

# Check if all files were processed by checking the error report
less 449532.stderr


#####################################


# Create bash script to merge stats info from .featureCounts from all samples
# into a single file:
for file in `find $HOME/IAV/count_summerisation/ \
-name *.featureCounts`; do echo echo \
"\`basename $file\` \`cut $file -f2 | sort | uniq -c | perl -p -e 's/\n/ /'\` >> \
annotation_summary_sense.txt" >> annotation_summary_sense.sh
done

nano annoation_summary_sense.sh
#!/bin/sh
#BSUB -J annotation_summary_job
#BSUB -n 10

bsub -q C < annoation_summary_sense.sh

# Check that all files were processed:
grep -c '.featureCounts' annotation_summary_sense.txt

# Copy all *sense-counts.txt files to temporary folder:
mkdir $HOME/IAV/count_summerisation/tmp

for file in `find $HOME/IAV/count_summerisation/ \
-name *sense-counts.txt`; do cp $file \
-t $HOME/IAV/count_summerisation/tmp; \
done

# Transfer all files from tmp to laptop, using FileZilla then remove tmp folder:
rm -r tmp

############################

have to rename


# Rename .featureCounts files:
for folder in \
`ls $HOME/IAV/Count_summarisation/Ion*`; \
do file2=`echo $folder`; \
echo "cd  $HOME/IAV/Count_summarisation/$file2; \
mv ./*_Aligned.out.sam.featureCounts ./${file2}_Aligned.sam.featureCounts" >> \
rename.sh; done;


for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/ \
-name *.featureCounts`; \
do outfile=`basename $file | perl -p -e \
's/home.ccorreia.scratch.PPDbRNAseqTimeCourse.STAR-2.5.1b_alignment.A*_*_*.//'`; \
mv $file $outfile >> rename.sh; done



########################################
# R analysis of gene counts with edgeR #
########################################

# Subsequent sense genes analyses were performed using the R statistical
# and the edgeR package. Please refer to file: PPDb-RNA-seq_paired_sense.R














