sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt-get update -qq
#sudo apt-get install -y cmake gcc gfortran g++ liblapack-dev libopenmpi-dev openmpi-bin
sudo apt-get install -y cmake gcc gfortran g++
sudo apt-get install -y netcdf-bin libnetcdf-dev
sudo apt-get install -y libopenmpi-dev openmpi-bin
sudo apt-get install -y libxml2 libxml2-dev libxml2-utils libxml-perl xml-core gnulib

which gcc
which g++
which gfortran
which ncdump
ls /usr/bin/mpi*

cd cime/scripts

./create_newcase --case f19_g16.ICLM45 --res f19_g16 --compset ICLM45 --mach travis-ci-linux --compiler gnu

cd f19_g16.ICLM45
./xmlchange DATM_CLMNCEP_YR_END=1972
./xmlchange PIO_TYPENAME=netcdf
./xmlchange RUNDIR=${PWD}/run
./xmlchange EXEROOT=${PWD}/bld
./xmlchange NTASKS=1
./xmlchange DIN_LOC_ROOT=$PWD
./case.setup
./case.build
cat bld/gptl.bldlog.*
