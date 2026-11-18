# BISICLES site specific notes

## Debian based systems (e.g Ubuntu and Ubuntu on WSL)

Debian GNU/Linux derivatives including Ubuntu 22.04 and 24.04 allow an
mpi environment , hdf5, and netcdf to be installed in a way that works
well for BISICLES. Likewise the Windows Subsystem for Linux (WSL), with
Ubuntu installed. On Ubuntu 22.04 or 24.04 run

      sudo apt install subversion build-essential g++ gfortran csh mpi-default-bin mpi-default-dev libhdf5-mpi-dev libhdf5-dev hdf5-tools libnetcdff-dev libnetcdf-dev python3 python3-dev libpython3-dev libfftw3-dev

If you want to use python for analysis you will need the usual libraries
(numpy, matplotlib etc) The file
 $BISICLES_HOME/Chombo/lib/mk/Make.defs.local should contain:

    PRECISION     = DOUBLE  
    CXX           = g++
    FC            = gfortran
    MPICXX        = mpiCC
    USE_HDF       = TRUE
    HDFINCFLAGS   = -I/usr/include/hdf5/serial/
    HDFLIBFLAGS   = -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5 -lz
    HDFMPIINCFLAGS= -I/usr/include/hdf5/openmpi/ 
    HDFMPILIBFLAGS= -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi/ -lhdf5  -lz
    cxxdbgflags    = -g -fPIC 
    cxxoptflags    = -fPIC -O2
    fdbgflags     =  -g -fPIC 
    foptflags     = -fPIC -O3 -ffast-math -funroll-loops
    USE_FFTW=TRUE
    FFTW_3=TRUE
    FFTWDIR=/usr

and the file  $BISICLES_HOME/BISICLES/code/mk/Make.defs. $mymachine
(where  $mymachine is the output from  'uname -n ') can be copied from
s $BISICLES_HOME/BISICLES/code/mk/Make.defs.ubuntu_22.4 (or 24.04) and
should contain

    PYTHON_INC=$(shell python3-config --includes)
    #--ldflags does not contain -lpython for reasons that escape me
    PYTHON_LIBS=-lpython3.10 $(shell python3-config --ldflags)
    NETCDF_HOME=$(shell nc-config --prefix)
    NETCDF_LIBS=-lnetcdff -lnetcdf -lhdf5_hl

or (for 24.04),

    PYTHON_INC=$(shell python3-config --includes)
    #--ldflags does not contain -lpython for reasons that escape me
    PYTHON_LIBS=-lpython3.12 $(shell python3-config --ldflags)
    NETCDF_HOME=$(shell nc-config --prefix)
    NETCDF_LIBS=-lnetcdff -lnetcdf -lhdf5_hl

## Rocky Linux 8

      sudo dnf install git subversion csh wget make gcc-c++ gcc-gfortran openmpi-devel python3-devel fftw3-devel perl zlib-devel
      module load mpi

Then follow the [generic instructions](readme.md). For BISICLES
configuration, use the file
 $BISICLES_HOME/BISICLES/code/mk/Make.defs.rocky8. For example,

      cd $BISICLES_HOME/BISICLES/code/mk
      ln -s Make.defs.rocky8 Make.defs.none

## University of Bristol

You might not need to compile BISICLES at UoB since a number of the
developers work there.

### Blue Pebble, GNU Compilers

Load the modules openmpi/5.0.3-et6p and languages/python/3.8.20 (the
more recent versions of python contain a libmpi version that conflicts
with the one in openmpi/5.0.3-et6p). Add the commands below to
 ~/.bash_profile.

      module load languages/python/3.8.20
      module load openmpi/5.0.3-et6p

Then follow the [generic instructions](readme.md). For BISICLES
configuration, used the file
 $BISICLES_HOME/BISICLES/code/mk/Make.defs.bp1. For example,

      cd $BISICLES_HOME/BISICLES/code/mk
      ln -s Make.defs.bp1 Make.defs.none

### Blue Crystal Phase 5, GNU Compilers

Blue Crystal Phase 5 in an x86 HPC managed by the Bristol Centre for Supercomputing. 
It is in early release, with [some documentation](https://docs.isambard.ac.uk/specs).
Some common programs are not yet installed, notably the is no python with 
development compinents (Python.h etc) and the system perl is missing
key packages for Chombo's build system. 

For now, hdf5, perl, python, and petsc can be found in /projects/b56av/bike, along
with Chombo and BISICLES. The bisicles binaries in /projects/b56av/bike/bisicles-uob
will be rebuilt often, but always from the master branch of https://github.com/CPOM-Land-Ice-Modelling/bisicles-uob.git.
Other bisicles-* directories are exprimental: use at your own risk.

### Isambard 3

The Isambard systems are UK HPC facilities managed by the Bristol Centre
for Supercomputing. See the [BriCS Isambard
documentation](https://docs.isambard.ac.uk) for details. BISICLES has
been built and runs on the main Isambard 3 Grace machine. This is an ARM
machine, built around NVIDIA Grace CPUs. BISICLES has been built and
runs well using the GNU compilers on this system. Isambard 3 Grace is a
Cray of some sort. To build and run, use the following modules.

    module load craype-network-ofi
    module load PrgEnv-Gnu
    module load cray-python/3.11.7

There is (as yet) no pre-installed hdf5 or PETSc, so follow the [generic
instructions](readme.md) to install those. There is no reason to build
serial hdf5. The file  $BISICLES_HOME/Chombo/lib/mk/Make.defs.local
should contain:

    MPI=TRUE 
    OPT=TRUE
    DEBUG=TRUE
    CXX=CC
    FC=ftn
    MPICXX=CC
    USE_64=TRUE
    CH_CPP=$(CXX) -E -P -C
    XTRACONFIG=.GNU
    cxxoptflags += -O3 -mcpu=neoverse-v2
    foptflags += -cpp -O3 -mcpu=neoverse-v2
    USE_HDF=TRUE
    # isambard 3 does not have an hdf5 module, so we need to build it.
    HDF5_PAR_DIR=...# path to parallel HDF5 install
    HDF5_SER_DIR=$(HDF5_PAR_DIR)# unless you decided to build serial versions for some reason 
    HDFLIBFLAGS=   -L$(HDF5_SER_DIR)/lib  -lhdf5 -lz -Wl,-rpath=$(HDF5_PAR_DIR)/lib
    HDFMPILIBFLAGS=-L$(HDF5_PAR_DIR)/lib  -lhdf5 -lz -Wl,-rpath=$(HDF5_PAR_DIR)/lib
    HDFINCFLAGS=   -I$(HDF5_SER_DIR)/include
    HDFMPIINCFLAGS=-I$(HDF5_PAR_DIR)/include

The file  $BISICLES_HOME/BISICLES/code/mk/Make.defs.none should include

    PYTHON_INC=$(shell python3-config --includes)
    #--ldflags does not include -lpython for reasons that escape me
    PYTHON_LIBS=-lpython3.11 $(shell python3-config --ldflags)
## NERSC

The US DOE-run National Energy Research Supercomputing Center (
[NERSC](http://www.nersc.gov)) has several machines which BISICLES users
may find useful.

### Cori

The Cray XC40 (
[Cori](https://www.nersc.gov/users/computational-systems/cori) ) at
NERSC offers Cray, GNU, and Intel compilers. Cori has two partitions  --
the Intel Xeon  "Haswell " partition is often faster execution
time-wise, but has longer waits in the queues. The Knights Landing
 "KNL " partition has slower run times, but often significantly faster
turnaround in the queues. I (Dan Martin) generally use the Haswell nodes
for things like visualization (I use VisIt in client-server mode,
leaving the data on NERSC), and the KNL nodes for moderate to large
production runs. We are also maintaining precompiled executables for our
own (and others ') use if that 's more convenient. If you have any
difficulties with these, e-mail Dan at DFMartin@lbl.gov.

#### Haswell nodes  -- precompiled binaries

The precompiled BISICLES executable for Haswell nodes (from the
public/trunk branch) can be found at

    /global/common/software/m1041/BISICLES/haswell/bin/driver2d.Linux.64.CC.ftn.OPT.MPI.PETSC.ex

To use it, you 'll need to load the python, GNU PrgEnv, and hdf5 modules
as below, and add the python library directory to your LD_LIBRARY_PATH
environment. In a csh/tcsh environment, this looks like:

    > module load craype-haswell
    > module swap PrgEnv-intel PrgEnv-gnu
    > module load cray-hdf5-parallel
    > module load python
    > module load cray-netcdf-hdf5parallel
    > module unload cray-shmem #needed for python only  
    > module load python
    > setenv LD_LIBRARY_PATH ${PYTHON_DIR}/lib:${LD_LIBRARY_PATH}

#### Haswell nodes  -- GNU compilers

To compile with GNU, switch to the correct PrgEnv and load an hdf5
module, plus netcdf and python if desired

    > module load craype-haswell
    > module swap PrgEnv-intel PrgEnv-gnu
    > module load cray-hdf5-parallel
    > module load python
    > module load cray-netcdf-hdf5parallel
    > module unload cray-shmem #needed for python only

To have the python interface work, the compiler needs the -shared and
-fPIC flags, and the linker needs the -dynamic flag. You may also need
to have the environment variable CRAYPE_LINK_TYPE set to  "dynamic "; in
the csh or tcsh shells:

    > setenv CRAYPE_LINK_TYPE dynamic

If you want to use the PETSc solvers, you can use the libraries
maintained by Mark Adams:

PETSC_DIR=/global/common/software/m1041/petsc_install/petsc_haswell_gnu

PETSC_ARCH= " "

You should be able to use the existing
Chombo/lib/mk/local/Make.defs.cori.hsw.gnu file in the Chombo release as
a starting point. An example
 $BISICLES_HOME/Chombo/lib/mk/local/Make.defs.local would be:

    makefiles+=Make.defs.local

    #default to MPI=TRUE,OPT=TRUE,DEBUG=FALSE 
    MPI=TRUE
    OPT=TRUE
    DEBUG=FALSE

    #this seems to be the Cray way
    CXX=CC
    FC=ftn
    MPICXX=CC
    USE_64=TRUE

    ifeq ($(PE_ENV),GNU)

    CH_CPP=$(CXX) -E -P 
    XTRACONFIG=.GNU
    cxxoptflags += -shared -fPIC 
    foptflags += -shared -fPIC
    ldoptflags += -dynamic
    cxxoptflags +=  -O3 -mavx2  -ffast-math 
    foptflags += -O3 -mavx2 -ffast-math 
    XTRALDFLAGS += -Wl,-zmuldefs

    else
    $(ECHO) "UNKNOWN PROGRAMMING ENVIRONMENT!"
    endif

    # The appropriate module (cray-hdf5-parallel) must be loaded for this to work.
    USE_HDF=TRUE
    HDFLIBFLAGS=   -L$(HDF5_DIR)/lib     $(HDF_POST_LINK_OPTS)  -lhdf5 -lz
    HDFMPILIBFLAGS=-L$(HDF5_DIR)/lib     $(HDF_POST_LINK_OPTS)  -lhdf5 -lz
    HDFINCFLAGS=   -I$(HDF5_DIR)/include $(HDF_INCLUDE_OPTS)
    HDFMPIINCFLAGS=-I$(HDF5_DIR)/include $(HDF_INCLUDE_OPTS)

For python and netcdf, the build should hopefully see the existing
Make.defs.cori file.

If you want to use the python functionality, then you should also load
the python module. Unfortunately, you currently also have to add the
python library directory to your LD_LIBRARY_PATH environment. In a
csh/tcsh environment, this looks like:

    module load python
    setenv LD_LIBRARY_PATH ${PYTHON_DIR}/lib:${LD_LIBRARY_PATH}

#### KNL Nodes  -- Intel compilers

The precompiled BISICLES executable for KNL nodes (from the public/trunk
branch) can be found at

    /global/common/software/m1041/BISICLES/KNL/bin/driver-public2d.Linux.64.CC.ftn.OPT.MPI.PETSC.ex

We 've found that the Intel compilers work well for the KNL part of
Cori. You 'll first need to load the correct modules and set some
environment variables (including pointing the build toward the Mark
Adams-maintained PETSC library).

    setenv NERSC_HOST `/usr/common/usg/bin/nersc_host`

    module unload PrgEnv-gnu
    module unload craype-haswell
    module load PrgEnv-intel
    module load craype-mic-knl
    module load cray-hdf5-parallel
    module load cray-netcdf-hdf5parallel
    module load visit  
    setenv PETSC_DIR PETSC_DIR=/global/common/software/m1041/petsc_install/petsc_knl_intel
    setenv PETSC_ARCH ""
    setenv CRAY_ROOTFS DSL
    setenv CRAYPE_LINK_TYPE dynamic

And then use Chombo/lib/mk/local/Make.defs.cori.knl If you want to use
the python functionality, then you should also load the python module.
Unfortunately, it appears that you currently also have to add the python
library directory to your LD_LIBRARY_PATH environment. In a csh/tcsh
environment, this looks like:

    module load python
    setenv LD_LIBRARY_PATH ${PYTHON_DIR}/lib:${LD_LIBRARY_PATH}

Your build should also automagically access the existing Make.defs.cori
file in BISICLES/code/mk, which is being kept up to date for KNL builds
and runs that we are constantly doing.

## ARCHER

The Cray XC30 at ARCHER (UK National HPC) offers Cray, GNU, and Intel
compilers. Cray is the default (but doesn 't seem to work for now). All
three compilers can be supported with a single Make.def.locals, it only
only necessary to load the correct PrgEnv module and the
cray-hdf5-parallel. One complexity is with statically versus dynamically
linked executables. Static linking is often simpler, but, will prevent
you from using the python interface as freely as you might like - in
particular,  'import math ', which is needed for many common math
functions, will result in a segfault. The alternative is dynamic
linking, which results in executables that have a longer list of
dependencies that must be matched between login nodes and compute nodes.
Choose static if you have no strong preference.

It is not generally a good idea to work with the filetools (nctoamr,
stats, flatten etc) on ARCHER - these are intended for interactive use,
and can be used on a workstation even with larger data, such as whole
Antarctic simulations. Ideally, you will have a GNU/linux workstation to
compile and run these tools, but if you have only a Windows or Mac OS X
machine, one possibility is to run a virtual machine with GNU/linux. We
have found that Oracle virtualbox and Ubuntu 18.04 work well on both
Windows and MacOS X machines.

### ARCHER  -- GNU compilers  -- static linking

For static linking, use the file Make.defs.archer included with Chombo

        ln -s $BISICLES_HOME/Chombo/lib/mk/local/Make.defs.archer $BISICLES_HOME/Chombo/lib/mk/Make.defs.local
      

To build statically linked executables with GNU, switch to the correct
PrgEnv and load an hdf5 module q

        > module unload PrgEnv-cray PrgEnv-intel
        > module load PrgEnv-gnuq
        > module load cray-hdf5-parallel
        > module load cray-petsc # if you want to build petsc-enabled executables
        > module load cray-netcdf-hdf5parallel # if you want to build the filetools
      

 $BISICLES_HOME/BISICLES/code/mk/Make.defs.archer includes the following
details for the python interface and netcdf.

       
        PYTHON_DIR=/work/y07/y07/cse/python/2.7.6-static/
        PYTHON_INC=-I$(PYTHON_DIR)/include/python2.7/
        PYTHON_LIBS=-L$(PYTHON_DIR)/lib -lpython2.7 -lpthread -ldl -lutil
        
        NETCDF_INC=-I$(NETCDF_DIR)/include
        NETCDF_LIBS=-L$(NETCDF_DIR)/lib -lnetcdf
      

Copy this file to Make.defs.none or Make.defs. $UNAMEN (where  $UNAMEN
is the output from uname -n) Then compile:

      cd BISICLES/code/exec2D
      make all OPT=TRUE MPI=TRUE DEBUG=FALSE # are all default 

the executable will be driver2d.Linux.64.CC.ftn.OPT.MPI.GNU.ex. To
compile a PETSc-enable executable,
driver2d.Linux.64.CC.ftn.OPT.MPI.PETSC.GNU.ex

         cd BISICLES/code/exec2D
         make all OPT=TRUE MPI=TRUE DEBUG=FALSE USE_PETSC=TRUE 
      

### ARCHER  -- GNU compilers  -- dynamic linking

For dynamic linking, use the file Make.defs.archer_dynamic_chombo
included with BISICLES

        ln -s $BISICLES_HOME/BISICLES/code/mk/Make.defs.archer_dynamic_chombo $BISICLES_HOME/Chombo/lib/mk/Make.defs.local
      

To build with GNU, switch to the correct PrgEnv and load an hdf5 module

        > module unload PrgEnv-cray PrgEnv-intel
        > module load PrgEnv-gnuq
        > module load cray-hdf5-parallel
        > module load python-compute
        > module load cray-netcdf-hdf5parallel
        > module load cray-petsc #if you want petsc
        > module unload cray-shmem
        > export CRAYPE_LINK_TYPE dynamic
      

 $BISICLES_HOME/BISICLES/code/mk/Make.defs.archer_dyamic includes the
following details for the python interface and netcdf.

        PYTHON_DIR=/work/y07/y07/cse/python/2.7.6/
        PYTHON_INC=-I$(PYTHON_DIR)/include/python2.7
        PYTHON_LIBS=-L$(PYTHON_DIR)/lib -lpython2.7 -lm -lpthread -lutil 

        NETCDF_INC=-I$(NETCDF_DIR)/include
        NETCDF_LIBS=-L$(NETCDF_DIR)/lib -lnetcdf

Copy this file to Make.defs.none or Make.defs. $UNAMEN (where  $UNAMEN
is the output from uname -n) Then compile

        cd BISICLES/code/exec2D
        make all OPT=TRUE MPI=TRUE DEBUG=FALSE # are all default 
      

the executable will be driver2d.Linux.64.CC.ftn.OPT.MPI.GNU.DY.ex. To
compile a PETSc-enable executable,
driver2d.Linux.64.CC.ftn.OPT.MPI.PETSC.GNU.DY.ex

        cd BISICLES/code/exec2D
        make all OPT=TRUE MPI=TRUE DEBUG=FALSE USE_PETSC=TRUE 
      

### ARCHER - Intel compilers

Compiling with Intel is essentially the same as with GNU. It 's just a
case of running

      > module unload PrgEnv-cray PrgEnv-gnu
      > module load PrgEnv-gnu
      > module load cray-hdf5-parallel # and so on

instead of

      > module unload PrgEnv-cray PrgEnv-intel
      > module load PrgEnv-gnu
      > module load cray-hdf5-parallel #and so on

The binary will be called driver2d.Linux.64.CC.ftn.OPT.MPI.INTEL.ex. (or
driver2d.Linux.64.CC.ftn.OPT.MPI.INTEL.DY.ex)

### Running jobs on ARCHER

ARCHER 's compute nodes do not have access to the /home mount point (and
therefore your home directory) It 's not obvious where login scripts
come from in that case, either. This can all be dealt with inside the
submission script (the MOM nodes do have access to home ...). A minimal
script (run from somewhere in /work) would look something like

    #PBS -l walltime=00:10:00 
    #PBS -j oe 
    #PBS -l select=1
    #PBS -A /your allocation/  

    export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)

    #load modules needed to get GLIBCXX, HDF5, python where the compute nodes can see them
    module unload PrgEnv-cray PrgEnv-intel
    module load PrgEnv-gnu
    module load cray-hdf5-parallel
    module load python-compute

    cd $PBS_O_WORKDIR

    #A GNU dynamic executable - change the name for intel, static linking, etc  
    EXE=/path/to/driver2d.Linux.64.CC.ftn.OPT.MPI.GNU.DY.ex 
    aprun -n 24 $EXE inputs.whatever

## The Monsoon Cray XC40 at the UK Met Office

The Monsoon Cray XC40 is much like ARCHER (and Edison and Cori at
NERSC). So much so that the ARCHER makefile seems to work, ie

    ln -s $BISICLES_HOME/Chombo/lib/mk/local/Make.defs.archer $BISICLES_HOME/Chombo/lib/mk/Make.defs.local

I 'm assuming that Monsoon users are interested in BISICLES coupled via
glimmer-cism to UKESM, i.e have checked out UniCiCles, and want to use
the intel compiler.

    > module swap PrgEnv-cray PrgEnv-intel
    > module load cray-hdf5-parallel
    > module load python/v2.7.9 #optional, needed if you want the python interface
    > module load cray-netcdf-hdf5parallel #optional, needed if you want glimmer-cism 
    > module load cray-tpsl/1.5.2 #optional, needed if you want the petsc solver
    > module load cray-petsc/3.6.1.0 #optional, needed if you want the petsc solver

If you want the python interface, make sure that
 $BISICLES_HOME/BISICLES/code/mk/Make.defs includes the following

          
    PYTHON_DIR=opt/python/gnu/2.7.9                                                                                    
    PYTHON_INC=$(PYTHON_DIR)/include/python2.7                                                                        
    PYTHON_LIBS=-L$(PYTHON_DIR)/lib -lpython2.7                                                            

and for netcdf stuff (e.g glimmer-cism, filetoools)

NETCDF_INC= $(NETCDF_DIR)/include NETCDF_LIBS=-L $(NETCDF_DIR)/lib
-lnetcdf

Compile with

    cd BISICLES/code/exec2D
    make all OPT=TRUE MPI=TRUE DEBUG=FALSE # are all default 

The binary will be called driver2d.Linux.64.CC.ftn.OPT.MPI.INTEL.ex

Or, with the PETSC interface

    cd BISICLES/code/exec2D
    make all OPT=TRUE MPI=TRUE DEBUG=FALSE USE_PETSC=TRUE

The binary will be called
driver2d.Linux.64.CC.ftn.OPT.MPI.PETSC.INTEL.ex
