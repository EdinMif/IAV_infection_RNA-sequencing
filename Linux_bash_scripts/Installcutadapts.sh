# Installation of Cutadapts1.12
# software can be downloaded via website https://pypi.python.org/pypi/cutadapt
#https://pypi.python.org/packages/41/9e/5b673f766dcf2dd787e0e6c9f08c4eea6f344ea8fce824241db93cc2175f/cutadapt-1.12.tar.gz

# download pip 
cd software
wget https://bootstrap.pypa.io/get-pip.py

# Before we can download pip we need to install Python-2.7 and override the existing 
# Python installed by the superuser to do that we need to do the following 
# Note that all the programs must be downloaded in the management node and not in the 

mkdir ~/python      
cd ~/python
wget https://www.python.org/ftp/python/2.7.11/Python-2.7.11.tgz
tar zxfv Python-2.7.11.tgz

        
#Install Python
./configure --prefix=$HOME/python
make 
make install 

mv python software/

# Create Symbolic link 
ln -s /home/ejmifsud/software/python/bin/python2.7 /home/ejmifsud/bin/python2.7

#call python
python2.7 --help

# download pip
#Perform the following tasks in management node

cd software
wget https://bootstrap.pypa.io/get-pip.py
python2.7 get-pip.py
rm get-pip.py

# Add to Path
ln -s /home/ejmifsud/software/python/bin/pip /home/ejmifsud/bin/pip

# Install cut adapt using pip
pip install --user --upgrade cutadapt 

#file is then downloaded and saved in following pathway
~/.local/bin/cutadapt --help
