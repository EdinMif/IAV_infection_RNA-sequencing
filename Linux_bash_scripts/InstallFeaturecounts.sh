# Installation of feature counts 



# Required package is featureCounts, which is part of Subread 1.5.0-p1 software,
# consult manual for details:
# http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf

# Download and install Subread 1.5.0-p1:
cd $HOME/Software
wget https://sourceforge.net/projects/subread/files/subread-1.5.1/subread-1.5.1-Linux-x86_64.tar.gz
tar zxvf subread-1.5.1-Linux-x86_64.tar.gz

# Remove compressed file
rm subread-1.5.1-Linux-x86_64.tar.gz

# create symbolic link
cd
ln -s /home/ejmifsud/software/subread-1.5.1-Linux-x86_64/bin/featureCounts /home/ejmifsud/bin/featureCounts
cd 
featureCounts --help
# if the program information appears it has been successfully installed

# NB: if the nano bash profile has not been set to bin then the following needs to be
 completed
nano .bash_profile
# Add :/home/ejmifsud/bin to the PATH=PATH:$HOME/bin
# Save it Press CNTL + X, then YES then enter. and restart session.
# To check the BASH_profile without editing type 
less .bash_profile 
# to exit press q