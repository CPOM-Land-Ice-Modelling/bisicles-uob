cd $BISICLES_HOME
echo `pwd`

#get hdf5 sources
if !(test -e hdf5-1.10.10.tar.bz2) then
    echo "downloading hdf5-1.10.9.tar.gz"
    wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.10/src/hdf5-1.10.10.tar.bz2
fi

mkdir -p hdf5/parallel/src
tar -jxf  hdf5-1.10.10.tar.bz2 -C hdf5/parallel/src

mkdir -p hdf5/serial/src
tar -jxf  hdf5-1.10.10.tar.bz2 -C hdf5/serial/src


#get netcdf sources

if !(test -e netcdf-c-4.9.2.tar.gz) then
    echo "downloading netcdf-4.9.2.tar.gz"
    wget https://downloads.unidata.ucar.edu/netcdf-c/4.9.2/netcdf-c-4.9.2.tar.gz
fi
mkdir -p netcdf/parallel/src
tar -zxf netcdf-c-4.9.2.tar.gz -C netcdf/parallel/src

mkdir -p netcdf/serial/src
tar -zxf netcdf-c-4.9.2.tar.gz -C netcdf/serial/src






