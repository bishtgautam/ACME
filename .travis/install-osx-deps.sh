# Install required software
brew update
brew install git
brew upgrade cmake
brew tap homebrew/science
brew unlink gcc
brew install gcc
brew install open-mpi
brew install netcdf
brew install make

# Make sure the weird gfortran library links are in place.
ln -s /usr/local/lib/gcc/5/libgfortran.dylib /usr/local/lib/libgfortran.dylib
ln -s /usr/local/lib/gcc/5/libgfortran.a /usr/local/lib/libgfortran.a

which gcc
which g++
which ncdump
gcc --version

ls -l /usr/local/bin/mpi*

cd cime/scripts
./create_newcase --case f19_g16.ICLM45 --res f19_g16 --compset ICLM45 --mach travis-ci-osx --compiler gnu

cd f19_g16.ICLM45
./xmlchange DATM_CLMNCEP_YR_END=1972
./xmlchange PIO_TYPENAME=netcdf
./xmlchange RUNDIR=${PWD}/run
./xmlchange EXEROOT=${PWD}/bld
./xmlchange NTASKS=1
./xmlchange DIN_LOC_ROOT=$PWD
./case.setup
alias gmake=make
cat Macros.make
cat env_mach_specific.xml
./case.build
cat bld/gptl.bldlog.*

