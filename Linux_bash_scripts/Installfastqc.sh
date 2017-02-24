## Install FastQC

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Download and install FastQC:
mkdir /home/ejmifsud/Software
# can not install fastqc from wget have to transfer file using filezilla to supercomputer
#transfered un-zipped file therefore do not need unzip command 

cd $HOME/software/FastQC
ls -lh
# no need to use chmod 755 fastqc function premissions already excutable 
/home/ejmifsud/software/FastQC/fastqc
 /home/ejmifsud/software/FastQC/fastqc --help


# Add it to PATH
mkdir /home/ejmifsud/bin
cd /home/ejmifsud/bin
ln -s /home/ejmifsud/software/FastQC/fastqc /home/ejmifsud/bin/fastqc
cd
nano .bash_profile
# Add :/home/ejmifsud/bin to the PATH=PATH:$HOME/bin
# Save it Press CNTL + X, then YES then enter. and restart session.
# To check the BASH_profile without editing type 
less .bash_profile 
# to exit press q



# Check if it worked:
cd 


fastqc --help