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

##############################################################################
# Alignment of FASTQ files against the Bos taurus reference genome with STAR #
##############################################################################

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


# Download annotation file for Mus_musculus genome version GRC38 Annotation Release-87 :
mkdir $HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/annotation_file

# exit internal node 
exit 
cd $HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/annotation_file
wget ftp://ftp.ensembl.org/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf.gz
 
# Log back into the Internal node 
 cd cd $HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/annotation_file
 gunzip Mus_musculus.GRCm38.87.gtf.gz
 
 
# Generate genome indexes files using annotations:
mkdir $HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/STAR-2.5.1b_index
cd $HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/STAR-2.5.1b_index

STAR --runThreadN 20 --runMode genomeGenerate \
--genomeDir /$HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/STAR-2.5.1b_index \
-- genomeFastaFiles \
$HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/source_file/Mus_musculus.GRCm38.dna.toplevel.fa \
--sjdbGTFfile /$HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/annotation_file/Mus_musculus.GRCm38.87.gtf \
--sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 99 \
--outFileNamePrefix  \
$HOME/IAV/genomes/mus_musculus/ensemble_GRCm38/STAR-2.5.1b_index




nohup STAR --runThreadN 20 --runMode genomeGenerate \
--genomeDir /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.1b_index \
--genomeFastaFiles \
/workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/source_file/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna \
--sjdbGTFfile /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/annotation_file/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff \
--sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 99 \
--outFileNamePrefix \
/workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.1b_index/Btau-UMD3.1.1 &

# Create and enter alignment working directory:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/STAR-2.5.1b_alignment
cd $HOME/scratch/PPDbRNAseqTimeCourse/STAR-2.5.1b_alignment

# Mapping reads from one FASTQ file to the indexed genome,
# to check if it works well:
nohup STAR --runMode alignReads --runThreadN 1 --genomeLoad LoadAndRemove \
--genomeDir /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.1b_index/ \
--readFilesIn \
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_001/trimmed_A6511_W10_F_R1_001.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_002/trimmed_A6511_W10_F_R1_002.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_003/trimmed_A6511_W10_F_R1_003.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_004/trimmed_A6511_W10_F_R1_004.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_005/trimmed_A6511_W10_F_R1_005.fastq.gz \
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_001/trimmed_A6511_W10_F_R2_001.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_002/trimmed_A6511_W10_F_R2_002.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_003/trimmed_A6511_W10_F_R2_003.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_004/trimmed_A6511_W10_F_R2_004.fastq.gz,\
$HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence/A6511_W10_F_005/trimmed_A6511_W10_F_R2_005.fastq.gz \
--readFilesCommand gunzip -c --outFilterMultimapNmax 10 \
--outFilterMismatchNmax 10 --outFileNamePrefix ./A6511_W10_F_ \
--outSAMtype BAM Unsorted --outReadsUnmapped Fastx &

# Create a bash script to perform alignment of paired FASTQ files:
for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/fastq_sequence \
-name *_R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/(_00.)/_002/g'`; \
file3=`echo $file | perl -p -e 's/(_00.)/_003/g'`; \
file4=`echo $file | perl -p -e 's/(_00.)/_004/g'`; \
file5=`echo $file | perl -p -e 's/(_00.)/_005/g'`; \
read1=`echo $file | perl -p -e 's/(R1_00.)/R2_001/'`; \
read2=`echo $file2 | perl -p -e 's/(R1_00.)/R2_002/'`; \
read3=`echo $file3 | perl -p -e 's/(R1_00.)/R2_003/'`; \
read4=`echo $file4 | perl -p -e 's/(R1_00.)/R2_004/'`; \
read5=`echo $file5 | perl -p -e 's/(R1_00.)/R2_005/'`; \
sample=`basename $file | perl -p -e 's/\_R1_001\.fastq\.gz//'`; \
foldername=`basename $sample | perl -p -e 's/trimmed\_//'`; \
echo "mkdir $HOME/scratch/PPDbRNAseqTimeCourse/STAR-2.5.1b_alignment/$foldername; \
cd $HOME/scratch/PPDbRNAseqTimeCourse/STAR-2.5.1b_alignment/$foldername; \
STAR --runMode alignReads --runThreadN 1 --genomeLoad LoadAndRemove \
--genomeDir /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.1b_index/ \
--readFilesIn $file,$file2,$file3,$file4,$file5 \
$read1,$read2,$read3,$read4,$read5 --readFilesCommand gunzip -c \
--outFilterMultimapNmax 10 --outFilterMismatchNmax 10 \
--outFileNamePrefix ./${foldername}_ --outSAMtype BAM Unsorted \
--outSAMattrIHstart 0 --outSAMattributes Standard --outReadsUnmapped Fastx" \
>> alignment.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 70 alignment.sh alignment.sh.
for script in `ls alignment.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check nohup.out file to see how many jobs finished successfully:
grep -c 'Finished successfully' alignment.sh.00.nohup
grep -c 'Finished successfully' alignment.sh.01.nohup

# Merge all STAR log.final.out files into a single file:
for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/STAR-2.5.1b_alignment \
-name *Log.final.out`; \
do perl /home/nnalpas/SVN/star_report_opener.pl -report $file; done;

#############################################
# FastQC quality check of aligned BAM files #
#############################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and go to working directory:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post_alignment
cd $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post_alignment

# Create a bash script to perform FastQC quality check on aligned SAM files:
for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/STAR-2.5.1b_alignment \
-name *.bam`; do echo "fastqc --noextract --nogroup -t 2 \
-o $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post_alignment $file" >> \
fastqc_aligned.sh; done;

# Split and run all scripts on Stampede
split -d -l 35 fastqc_aligned.sh fastqc_aligned.sh.
for script in `ls fastqc_aligned.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Delete all the HTML files:
rm -r *.html

# Check all output from FastQC:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post_alignment/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post_alignment/tmp; \
done

for file in \
`find $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post_alignment/tmp \
-name summary.txt`; do more $file >> reports_post-alignment.txt; done

for file in \
`find $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/post_alignment/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_post_alignment.txt; \
done

# Check if all files were processed:
grep -c '##FastQC' basic_stats_post_alignment.txt
grep -c 'Basic Statistics' reports_post-alignment.txt
grep -c 'Analysis complete' fastqc_aligned.sh.00.nohup
grep -c 'Analysis complete' fastqc_aligned.sh.01.nohup
grep -c 'Analysis complete' fastqc_aligned.sh.02.nohup
grep -c 'Analysis complete' fastqc_aligned.sh.03.nohup

# Remove temporary folder:
rm -r tmp/

###################################################################
# Summarisation of gene counts with featureCounts for sense genes #
###################################################################

# Required package is featureCounts, which is part of Subread 1.5.0-p1 software,
# consult manual for details:
# http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf

# Create working directories:
cd $HOME/scratch/PPDbRNAseqTimeCourse/
mkdir -p Count_summarisation/sense
cd $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense

# Run featureCounts with one sample to check if it is working fine:
featureCounts -a \
/workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/annotation_file/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff \
-B -p -C -R -s 1 -T 15 -t gene -g Dbxref -o ./counts.txt \
$HOME/scratch/PPDbRNAseqTimeCourse/STAR-2.5.1b_alignment/A6511_W10_F/A6511_W10_F_Aligned.out.bam

# Create a bash script to run featureCounts on BAM file containing multihits and
# uniquely mapped reads using the reversely stranded parameter:
for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/STAR-2.5.1b_alignment \
-name *_Aligned.out.bam`; \
do sample=`basename $file | perl -p -e 's/_Aligned.out.bam//'`; \
echo "mkdir $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/$sample; \
cd $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/$sample; \
featureCounts -a \
/workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/annotation_file/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff \
-B -p -C -R -s 1 -T 10 -t gene -g Dbxref \
-o ${sample}_sense-counts.txt $file" >> sense_count.sh; done

# Split and run all scripts on Stampede:
split -d -l 70 sense_count.sh sense_count.sh.
for script in `ls sense_count.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all files were processed:
grep -c 'Read assignment finished.' sense_count.sh.00.nohup
grep -c 'Read assignment finished.' sense_count.sh.01.nohup


############################

have to rename


# Rename .featureCounts files:
for folder in \
`ls $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/A6*`; \
do file2=`echo $folder`; \
echo "cd $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/$file2; \
mv ./*_Aligned.out.sam.featureCounts ./${file2}_Aligned.sam.featureCounts" >> \
rename.sh; done;


for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/ \
-name *.featureCounts`; \
do outfile=`basename $file | perl -p -e \
's/home.ccorreia.scratch.PPDbRNAseqTimeCourse.STAR-2.5.1b_alignment.A*_*_*.//'`; \
mv $file $outfile >> rename.sh; done


#####################################


# Create bash script to merge stats info from .featureCounts from all samples
# into a single file:
for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/ \
-name *.featureCounts`; do echo echo \
"\`basename $file\` \`cut $file -f2 | sort | uniq -c | perl -p -e 's/\n/ /'\` >> \
annotation_summary_sense.txt" >> annotation_summary_sense.sh
done

# Split and run scripts on Stampede:
split -d -l 70 annotation_summary_sense.sh annotation_summary_sense.sh.
for script in `ls annotation_summary_sense.sh.*`
do
chmod 755 $script
nohup ./$script &
done

# Check that all files were processed:
grep -c '.featureCounts' annotation_summary_sense.txt

# Copy all *sense-counts.txt files to temporary folder:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/tmp

for file in `find $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/ \
-name *sense-counts.txt`; do cp $file \
-t $HOME/scratch/PPDbRNAseqTimeCourse/Count_summarisation/sense/tmp; \
done

# Transfer all files from tmp to laptop, using WinSCP, then remove tmp folder:
rm -r tmp


########################################
# R analysis of gene counts with edgeR #
########################################

# Subsequent sense genes analyses were performed using the R statistical
# and the edgeR package. Please refer to file: PPDb-RNA-seq_paired_sense.R














