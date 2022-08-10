# SPARCED (non official) installation guide on Ubuntu
---
Hi! 🌄

If you're new to SPARCED and wish to get a working environment setup on Ubuntu, then you're at the right place!
This is a document I wrote as a summer intern at the Birtwistle lab.
For some reason my previous setup was buggy and I had to reinstall everything from the ground up,
so I decided to write this short document to make the process easier for newcomers like me 🙂

## Environment
I am running an Ubuntu 22.04 LTS VM on VirtualBox.

Run the following command to make sure everything is up to date:
```bash
sudo apt-get update
sudo apt-get upgrade
```

## VirtualBox Guest Additions
_If you are not using VirtualBox, skip this section._

Insert the guest additions CD, then run:
```bash
sudo apt install build-essentials dkms linux-headers-generic # Necessary for minimal installations of Ubuntu
lsblk | grep "rom"
cd /media/{user}/VBox_GAs_6.1.34 # with {user} being your username
# Please note that the version number might differ but this is usually not an issue
sudo ./VBoxLinuxAdditions.run
```
You will need to restart the VM after that. Don't forget to eject the CD! 😉

## Git, GitHub & SSH
```bash
# Git
sudo apt install git-all
# SSH for GitHub
ssh-keygen -t ed25519 -C "{email}" # with {email} being your email for GitHub
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519
cat ~/.ssh/id_ed25519 # Prints your SSH key in the terminal
```
Copy your SSH key and in your GitHub settings, create a new SSH key where you can paste it.

Run the following command to test your SSH connexion:
```bash
ssh -T git@github.com
```

## Anaconda
Download the [Anaconda installer for Linux](https://www.anaconda.com/products/distribution#linux), then run the following commands:
```bash
sudo apt-get install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6
bash ~/Downloads/Anaconda3-{version-number}-Linux-x86_64.sh # with {version-number} being your version number
source ~/.bashrc # or alternatively close and reopen your terminal
conda config --set auto_activate_base False # set it according to your preferences
```
Verify your installation using:
```bash
conda list
```

## Packages
```bash
conda create -n sparced # Creates an environment named "sparced"
source activate sparced # Activates the "sparced" environment
conda install -c anaconda spyder
conda install matplotlib
conda install pandas
pip install -Iv antimony==2.12.0.1 # WARNING: antimony >= 2.13.0 doesn't work with SPARCED
pip install python-libsbml
pip install scipy
```
### The Amici Package
```bash
sudo apt install libatlas-base-dev swig
pip install amici
```
You might get an error about the CBLAS library (this happens mostly on Palmetto), to fix it run:
```bash
conda install -c conda-forge openblas
export BLAS_LIBS=-lopenblas
```
### TODO: mpi4py
_If you are not going to use parallel computation, skip this section._

## TODO: Docker
_If you are not going to run SPARCED inside the official Jupyter Notebook container, skip this section._

## SPARCED 🎆
This is only a setup suggestion:
```bash
cd Documents
mkdir birtwistle-lab ; cd birtwistle-lab
git clone --recursive https://github.com/birtwistlelab/SPARCED.git # The official SPARCED repository
cd ..
git clone --recursive https://github.com/{username}/SPARCED.git # with {username} being your username on GitHub, assuming that you already forked SPARCED
git clone --recursive https://github.com/ChocolateCharlie/sparced-personal-utilities.git # My code, feel free to improve it :)
```

## OpenMPI
_If you are not going to use parallel computation, skip this section._

You will need OpenMPI if you want to do some parallel computing (or mostly debug on your own machine some code intended to run parallely on Palmetto).

First, download the latest stable version of OpenMPI from the [openmpi.org](https://www.open-mpi.org//software/ompi/v4.1/) website. You want the ```.tar.gz``` extension.
Then run the following commands (some can take a few minutes and be very verbose, so get a ☕):
```bash
mkdir openmpi
cd openmpi
cp ~/Downloads/openmpi-1.8.7.tar.gz # or whatever version number you have
tar -xzvf openmpi-1.8.7.tar.gz # or whaterver version number you have
cd openmpi-1.8.7 # or whaterver version number you have
./configure --prefix=$HOME/openmpi
make install
export PTH=$HOME/openmpi/bin:$PATH
export LD_LIBRARY_PATH=$HOME/openmpi/lib:$LD_LIBRARY_PATH
```
You can check if the installation process worked using:
```bash
mpirun --version
```
:warning: If you get an error involving Fortran during this process and you find a way to fix it without having to reinstall everything, please notify me!

## Clean
Remove all unused packages that were installed by dependencies during the setup:
```bash
sudo apt-get autoremove
```

Congratulations! You now have a full setup of SPARCED! 🦠
