# Required software is STAR 2.5.2b, consult manual/tutorial for details:
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

# Download and install STAR:
https://github.com/alexdobin/STAR/archive/2.5.2b.tar.gz

#Downloaded file was transfered to CZC supercomputer using Filezilla

# Add it to PATH
cd
ln -s /home/ejmifsud/software/STAR-2.5.2b/bin/Linux_x86_64/STAR /home/ejmifsud/bin/STAR
cd

#Call the program 
STAR --help
# Once the program information comes up it is correctly installed 
# restart session.


# Install STARlong used with PACBIO sequences
# program was already downloaded as part of the Linux_x86_64 package therefore we only
# need to create a symbolic link to the bin folder

# Creating a Symbolic link
cd
ln -s /home/ejmifsud/software/STAR-2.5.2b/bin/Linux_x86_64/STARlong \
/home/ejmifsud/bin/STARlong
cd

STARlong --help


# NB: if the nano bash profile has not been set to bin then the following needs to be
 completed
nano .bash_profile
# Add :/home/ejmifsud/bin to the PATH=PATH:$HOME/bin
# Save it Press CNTL + X, then YES then enter. and restart session.
# To check the BASH_profile without editing type 
less .bash_profile 
# to exit press q

# Ignore all this (pathway to the programes) 
home/ejmifsud/software/STAR-2.5.2b/bin/Linux_x86_64
home/ejmifsud/software/STAR-2.5.2b/bin/Linux_x86_64/STAR
home/ejmifsud/software/STAR-2.5.2b/bin/Linux_x86_64/STARlong
