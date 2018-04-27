# Install required software
brew update
brew cask uninstall --force oclint
brew upgrade cmake
brew install gcc@7
brew install mpich
brew install netcdf --with-fortran
brew install make

# Make sure the weird gfortran library links are in place.
ln -s /usr/local/lib/gcc/5/libgfortran.dylib /usr/local/lib/libgfortran.dylib
ln -s /usr/local/lib/gcc/5/libgfortran.a /usr/local/lib/libgfortran.a

which gcc
which g++
which ncdump
gcc --version

ls -l /usr/local/bin/mpi*
ls -l /usr/local/lib/libnetcdf*

/usr/local/bin/mpif90 --version

which nc-config

nc-config --all

nc-config --flibs

