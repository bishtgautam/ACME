#!/bin/sh

export PETSC_DIR=$PWD/petsc
export PETSC_ARCH=petsc-arch


cd $PETSC_DIR/$PETSC_ARCH/externalpackages

wget http://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-fortran-4.4.1.tar.gz

tar -xzvf netcdf-fortran-4.4.1.tar.gz
cd netcdf-fortran-4.4.1
export FC=$PETSC_DIR/$PETSC_ARCH/bin/mpif90 
export F77=$PETSC_DIR/$PETSC_ARCH/bin/mpif77
export CC=$PETSC_DIR/$PETSC_ARCH/bin/mpicc

export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
export CPPFLAGS="-DgFortran -I$PETSC_DIR/$PETSC_ARCH/include"

export LDFLAGS=-L$PETSC_DIR/$PETSC_ARCH/lib
./configure --prefix=$PETSC_DIR/$PETSC_ARCH

make 
make check
make install


cd $PETSC_DIR/$PETSC_ARCH/bin
./nc-config --all

ls -l $PETSC_DIR/$PETSC_ARCH/include

ls -l $PETSC_DIR/$PETSC_ARCH/lib/lib*

ls -l $PETSC_DIR/$PETSC_ARCH/lib/libnetcdff.a

which nm

nm $PETSC_DIR/$PETSC_ARCH/lib/libnetcdff.a

nm $PETSC_DIR/$PETSC_ARCH/lib/libnetcdff.a | grep deflate