sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt-get update -qq
sudo apt-get install -y cmake gcc-5 gfortran-5 g++-5 liblapack-dev libopenmpi-dev
#sudo apt-get install -y libnetcdf-dev netcdf-bin
sudo apt-get install -y libxml2 libxml2-dev libxml2-utils libxml-perl xml-core gnulib

#sudo apt-get install -y libmpich-dev mpich

#sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
#sudo apt-get update -qq
#sudo apt-get install -y cmake gcc gfortran g++ liblapack-dev libopenmpi-dev openmpi-bin
#sudo apt-get install -y cmake gcc-5 gfortran-5 g++-5
#sudo apt-get install -y libopenmpi-dev openmpi-bin
#sudo apt-get install -y cmake liblapack-dev
#sudo apt-get install -y netcdf-bin libnetcdf-dev
#sudo apt-get install -y libmpich-dev mpich
#sudo apt-get install -y libxml2 libxml2-dev libxml2-utils libxml-perl xml-core gnulib




which gcc
which g++
which gfortran
which ncdump
ls -l /usr/bin/mpi*
ls -l /usr/bin/gcc*

.travis/install-petsc.sh

ls -l ${PWD}/petsc/petsc-arch/bin/
${PWD}/petsc/petsc-arch/bin/nc-config --all
ls -l ${PWD}/petsc/petsc-arch/externalpackages
head -20 ${PWD}/petsc/petsc-arch/externalpackages/mpich-3.2/config.log
head -20 ${PWD}/petsc/petsc-arch/externalpackages/netcdf-4.4.1/config.log
.travis/install-netcdff.sh
cd cime/scripts
./create_newcase --case f19_g16.ICLM45 --res f19_g16 --compset ICLM45 --mach travis-ci-$TRAVIS_OS_NAME --compiler gnu
cd f19_g16.ICLM45
./xmlchange DATM_CLMNCEP_YR_END=1972
./xmlchange PIO_TYPENAME=netcdf
./xmlchange RUNDIR=${PWD}/run
./xmlchange EXEROOT=${PWD}/bld
./xmlchange NTASKS=1
./xmlchange DIN_LOC_ROOT=$PWD
./xmlchange MPILIB=mpich
./case.setup
cat Macros.make
./case.build
cat bld/gptl.bldlog.*
cat bld/pio.bldlog.*
cat bld/lnd.bldlog*


#cd cime/scripts
#
#./create_newcase --case f19_g16.ICLM45 --res f19_g16 --compset ICLM45 --mach travis-ci-linux --compiler gnu
#
#cd f19_g16.ICLM45
#./xmlchange DATM_CLMNCEP_YR_END=1972
#./xmlchange PIO_TYPENAME=netcdf
#./xmlchange RUNDIR=${PWD}/run
#./xmlchange EXEROOT=${PWD}/bld
#./xmlchange NTASKS=1
#./xmlchange DIN_LOC_ROOT=$PWD
#./case.setup
#cat Macros.make
#./case.build
#cat bld/gptl.bldlog.*
#cat bld/lnd.bldlog*
