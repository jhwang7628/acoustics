#!/bin/bash

cd $HOME
mkdir code
sudo mkdir /hdd1/$USER
sudo chown -R $USER /hdd1/$USER
mkdir /hdd1/$USER/data
ln -s /hdd1/$USER/data $HOME/data
sudo apt-get install git -y

cd code
git config --global user.name "Ante Qu"
git config --global user.email antequ@cs.stanford.edu
git config --global core.editor vim
git clone https://antequ@bitbucket.org/jw969/acoustics.git

cd acoustics
git fetch && git checkout dev_antequ 
git submodule init
git submodule update

#install krb5-user krb5-doc
#to install:
# MKL
# boost
# Eigen (after superlu gmp glut?)
# CGAL
# gmp
# libconfig
# SuperLU
# VLFeat
# libqglviewer
# gts?
# libsndfile?
# glut?
#gts

sudo apt-get install libqt4-dev qt4-doc -y
sudo apt-get install cmake qt5-default gcc python2.7 paraview -y
sudo apt-get install qt5-doc gcc-doc python2.7-doc cmake-doc -y
sudo apt-get install nvidia-cuda-dev nvidia-cuda-doc nvidia-cuda-toolkit -y
# MKL
# use student license
sudo apt-get install lib32stdc++6 lib32stdc++-5-dev lib32gcc-5-dev gcc-multilib -y
read -p "Please install Intel MKL. Press enter to continue"
#BOOST
sudo apt-get install libboost-all-dev -y
#GMP and mpfr
sudo apt-get install libgmp-dev libmpfr-dev libmpfr-doc gmp-doc libgmp10-doc -y
#libsuperlu
sudo apt-get install libsuperlu-dev libsuperlu-doc -y
sudo apt-get install gfortran-doc gfortran-5-doc liblapack-doc-man liblapack-doc -y
#FFTW
sudo apt-get install libfftw3-dev libfftw3-doc -y
#GLUT
sudo apt-get install freeglut3-dev -y
#libqglviewer
sudo apt-get install libqglviewer-dev libqglviewer-doc -y
#libconfig
sudo apt-get install libconfig-dev libconfig-doc libconfig++-dev -y
#cgal (should I get libcgal-dev?)
sudo apt-get install libcgal-qt5-dev -y
#gts
sudo apt-get install libgts-dev -y
#vlfeat
sudo apt-get install libvlfeat-dev libvlfeat-doc -y
#libsndfile
sudo apt-get install libsndfile1-dev -y
#gsl
sudo apt-get install libgsl-dev gsl-doc-info gsl-doc-pdf gsl-ref-psdoc gsl-ref-html -y
#arpack
sudo apt-get install libarpack2-dev libarpack++2-dev -y
#xi
sudo apt-get install libxi-dev libxi6 -y
#iomp5
sudo apt-get install libiomp-dev libiomp-doc -y
mkdir $HOME/installers
pushd $HOME/installers
# EIGEN INSTALLATION
sudo apt-get install libeigen3-dev libeigen3-doc -y
#echo "Installing Eigen..."
#wget http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz
#mv 3.3.4.tar.gz eigen3.3.4.tar.gz
#tar xvzf eigen3.3.4.tar.gz
#cd eigen-eigen-5a0156e40feb
#mkdir build
#cd build
#cmake ..

#make acoustics
pushd

echo "export EIGEN_ROOT=/usr/include/eigen3" >> ~/.bashrc
export EIGEN_ROOT=/usr/include/eigen3
mkdir build_release
cd build_release
cmake ..

make -j36 rigidsim-gui rigidsim
make -j36 unit_testing_FDTD_AcousticSimulator
make -j36 fdtd-acoustic-simulator-viewer
make -j36 tetviewer
make -j36 isostuffer
echo "export PATH=\$PATH:\$HOME/code/acoustics/build_release/bin" >> ~/.bashrc
export PATH=$PATH:$HOME/code/acoustics/build_release/bin

echo "export SCRIPTS=\$HOME/code/acoustics/src/scripts" >> ~/.bashrc
export SCRIPTS=$HOME/code/acoustics/src/scripts
echo "export PATH=\$PATH:\$SCRIPTS" >> ~/.bashrc
export PATH=$PATH:$SCRIPTS

#Mitsuba
#https://www.mitsuba-renderer.org/releases/current/documentation.pdf
read -p 'Press [Enter] key to continue to Mitsuba installation...'
pushd $HOME/installers
sudo apt-get install build-essential scons mercurial qt4-dev-tools libpng12-dev libjpeg-dev libilmbase-dev libxerces-c-dev libboost-all-dev libopenexr-dev libglewmx-dev libxxf86vm-dev libpcrecpp0v5 libeigen3-dev libfftw3-dev -y
mkdir mitsuba
cd mitsuba
# these don't work on 16.04
#wget http://www.mitsuba-renderer.org/releases/current/trusty/collada-dom-dev_2.4.0-1_amd64.deb
#wget http://www.mitsuba-renderer.org/releases/current/trusty/collada-dom_2.4.0-1_amd64.deb
#wget http://www.mitsuba-renderer.org/releases/current/trusty/mitsuba-dev_0.5.0-1_amd64.deb
#wget http://www.mitsuba-renderer.org/releases/current/trusty/mitsuba_0.5.0-1_amd64.deb

#sudo dpkg --install collada-dom_*.deb
#sudo dpkg --install mitsuba*.deb

#following these instead: 
# https://cgcvtutorials.wordpress.com/2017/05/31/install-mitsuba-on-linux-16-04/

wget https://www.mitsuba-renderer.org/repos/mitsuba/archive/tip.zip

sudo apt-get install unzip -y
unzip tip.zip
mkdir $HOME/installed

mv mitsuba-af602c6fd98a $HOME/installed/
cd $HOME/installed/mitsuba-af602c6fd98a/

cp build/config-linux-gcc.py config.py

sudo apt-get install scons qt4-dev-tools libboost-all-dev libglewmx-dev libpcrecpp0v5 libxerces-c-dev -y

echo "Please follow the instructions in https://cgcvtutorials.wordpress.com/2017/05/31/install-mitsuba-on-linux-16-04/"
echo "and fix src/bsdfs/irawan.h as well as include/mitsuba/core/autodiff.h."
read -p "Press [Enter] when done"
scons –j 36

echo "loadmitsuba() {" >> ~/.bashrc
echo "  . ~/installed/mitsuba-af602c6fd98a/setpath.sh" >> ~/.bashrc
echo "}" >> ~/.bashrc



