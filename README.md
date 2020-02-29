## OpenFoam 6.0 install on SuperMike II
Follow the guide @ http://openfoamwiki.net/index.php/Installation/Linux/OpenFOAM-6/CentOS_SL_RHEL#CentOS_6.10_.28without_SCL.29
```bash

#Using node for faster compile
qsub -I -l nodes=2:ppn=16,walltime=00:120:00 -q workq -A hpc_startup_xjwang

#enable module
sh /usr/local/packages/Modules/module_setup.sh

#Download and unpack
cd /project/xjwang/
mkdir OpenFOAM
cd OpenFOAM
git clone https://github.com/OpenFOAM/OpenFOAM-6.git
git clone https://github.com/OpenFOAM/ThirdParty-6.git

cd ThirdParty-6
mkdir download

wget -P download https://www.cmake.org/files/v3.9/cmake-3.9.0.tar.gz
wget -P download  https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.10/CGAL-4.10.tar.xz
wget -P download https://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.bz2
wget -P download https://www.open-mpi.org/software/ompi/v2.1/downloads/openmpi-2.1.1.tar.bz2

tar -xzf download/cmake-3.9.0.tar.gz
tar -xJf download/CGAL-4.10.tar.xz
tar -xjf download/boost_1_55_0.tar.bz2
tar -xjf download/openmpi-2.1.1.tar.bz2

cd ..

#Change Default software version

sed -i -e 's/\(boost_version=\)boost-system/\1boost_1_55_0/' OpenFOAM-6/etc/config.sh/CGAL
sed -i -e 's/\(cgal_version=\)cgal-system/\1CGAL-4.10/' OpenFOAM-6/etc/config.sh/CGAL

#Check arch version, it should be x86_64 in SuperMike II
uname -m

export PROJECT=/project/xjwang
printf '\n%s\n' 'export PROJECT=/project/xjwang' >> $HOME/.bashrc


#Setup building parameters
source $PROJECT/OpenFOAM/OpenFOAM-6/etc/bashrc WM_COMPILER_TYPE=ThirdParty WM_COMPILER=Gcc49 WM_LABEL_SIZE=64 WM_MPLIB=OPENMPI FOAMY_HEX_MESH=yes

#For faster compile using 16 cores in SuperMike II
export WM_NCOMPPROCS=16

#Then save an alias in the personal .bashrc file, simply by running the following command:
echo "alias of6='source \$PROJECT/OpenFOAM/OpenFOAM-6/etc/bashrc $FOAM_SETTINGS'" >> $HOME/.bashrc

#Build GCC and cmake
cd $WM_THIRD_PARTY_DIR
wget "https://raw.github.com/wyldckat/scripts4OpenFOAM3rdParty/master/getGcc"
wget "https://raw.github.com/wyldckat/ThirdParty-2.0.x/binutils/makeBinutils"
wget "https://raw.github.com/wyldckat/ThirdParty-2.0.x/binutils/getBinutils"
chmod +x get* make*


cd $WM_THIRD_PARTY_DIR
./makeCmake > log.makeCmake 2>&1
wmRefresh

./getGcc gcc-4.9.3 gmp-5.1.2 mpfr-3.1.2 mpc-1.0.1
./makeGcc -no-multilib > log.makeGcc 2>&1
wmRefresh

#build a custom GNU Binutils
./getBinutils
./makeBinutils gcc-4.9.3 > log.makeBinutils 2>&1

#Build CGAL

#Fix Gcc Error in case of gcc 4.8.5 is used 
#module load gcc/4.9.2
#cd $PROJECT/OpenFOAM/ThirdParty-6/platforms/linux64/gcc-4.8.5/lib64/
#mv libstdc++.so.6 libstdc++.so.6.orig
#ln -s /home/compilers/gcc/4.9.2/lib64/libstdc++.so.6 libstdc++.so.6

./makeCGAL > log.makeCGAL 2>&1
wmRefresh

#Build OpenMPI
cd $WM_THIRD_PARTY_DIR
./Allwmake -j $WM_NCOMPPROCS > log.make 2>&1
wmRefresh

#------!!!Build OpenFoam!!!--------
cd $WM_PROJECT_DIR
./Allwmake -j $WM_NCOMPPROCS > log.make 2>&1

#------!!!Testing!!!--------
of6
icoFoam -help


```
Reminder: Whenever you start a new terminal window or tab, you should run the alias command associated to the OpenFOAM 6 shell environment. In other words, run the following command whenever you start a new terminal:
```
of6
```
