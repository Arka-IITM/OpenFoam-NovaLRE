#!/bin/bash
#PBS -o logfile.log
#PBS -e errorfile.err
#PBS -l walltime=00:20:00
#PBS -l select=1:ncpus=40


tpdir=`echo $PBS_JOBID | cut -f 1 -d .`
tempdir=$HOME/scratch/job$tpdir
mkdir -p $tempdir
cd $tempdir
exec > "$tempdir/master.log" 2>&1
set -x  
export PS4='+ $(date +%T) '   # prefix each traced line with HH:MM:SS
echo "copying" >> log
mv $PBS_O_WORKDIR/* .          # copies system/ along with everything else


module load openfoam2406
source /lfs/sware/OpenFOAM/OpenFOAM-2406/etc/bashrc
export PATH=/lfs/sware/gcc13.1.0/bin:/lfs/sware/openmpi405/bin:$PATH
export LD_LIBRARY_PATH=/lfs/sware/gcc13.1.0/lib64:/lfs/sware/gcc13.1.0/lib:/lfs/sware/openmpi405/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=/lfs/sware/gcc13.1.0/lib64:/lfs/sware/gcc13.1.0/lib:$LIBRARY_PATH
export WM_CC=/lfs/sware/gcc13.1.0/bin/gcc
export WM_CXX=/lfs/sware/gcc13.1.0/bin/g++

wclean applications/atomizationFoam
wclean src/libAtomization
rm -f /lfs/usrhome/btech/ed25b061/OpenFOAM/ed25b061-v2406/platforms/linux64GccDPInt32Opt/lib/libAtomization.so

export WM_CFLAGS="-L/lfs/sware/gcc13.1.0/lib64 -Wl,-rpath,/lfs/sware/gcc13.1.0/lib64"
export WM_CXXFLAGS="-L/lfs/sware/gcc13.1.0/lib64 -Wl,-rpath,/lfs/sware/gcc13.1.0/lib64"
export WM_LDFLAGS="-L/lfs/sware/gcc13.1.0/lib64 -Wl,-rpath,/lfs/sware/gcc13.1.0/lib64"


# Strip gcc9.2.0 / Intel / cuda include pollution so gcc13 uses its OWN headers
unset CPATH CPLUS_INCLUDE_PATH C_INCLUDE_PATH GCC_EXEC_PREFIX

echo "Allwmake" >> log
which g++
g++ --version          # must say 13.1.0, NOT 9.2.0
echo "$WM_CC $WM_CXX"a > log.allwmake
export WM_CC=/lfs/sware/gcc13.1.0/bin/gcc
export WM_CXX=/lfs/sware/gcc13.1.0/bin/g++

which g++
g++ --version          # must say 13.1.0, NOT 9.2.0
echo "$WM_CC $WM_CXX" >> log.allwmake

./Allwmake >> log.allwmake
echo "ran allwmake" >> log

export PATH=/lfs/usrhome/btech/ed25b061/OpenFOAM/ed25b061-v2406/platforms/linux64GccDPInt32Opt/bin:$PATH
export LD_LIBRARY_PATH=/lfs/usrhome/btech/ed25b061/OpenFOAM/ed25b061-v2406/platforms/linux64GccDPInt32Opt/lib:$LD_LIBRARY_PATH


cd run || { echo "ERROR: run/crossFlow not found under $tempdir" >> log; exit 1; }
[ -f system/controlDict ]        || { echo "ERROR: system/controlDict missing in $(pwd)" >> log; exit 1; }
[ -f system/decomposeParDict ]   || { echo "ERROR: system/decomposeParDict missing in $(pwd)" >> log; exit 1; }

# ./Allrun
echo "Allrun complete" >> log.allruncomplete

#echo "decompose" >> log
 decomposePar > log.decomposePar 2>&1
#echo "done" >> log

echo "Starting atomisationFoam" >> log.startingFoam
mpirun --oversubscribe --use-hwthread-cpus -np 40 atomisationFoam -parallel -fileHandler masterUncollated > log.Foam 2>&1 || mpirun --oversubscribe --use-hwthread-cpus -np 40 \
  /lfs/usrhome/btech/ed25b061/OpenFOAM/ed25b061-v2406/platforms/linux64GccDPInt32Opt/bin/atomizationFoam \
  -parallel -fileHandler masterUncollated > log.Foam 2>&1
echo "Calculation done. Reconstructing" >> log.startingFoam
reconstructPar > log.reconstructPar 2>&1
