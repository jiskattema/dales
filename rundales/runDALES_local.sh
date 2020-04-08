#!/bin/bash
#----------------------------------------------------------------------
# Script for running DALES#
#----------------------------------------------------------------------
#
# - requires following variables:
#     $USERN    ::  username
#     $DAYC      ::  current day of model of code adjustment/compilation
#  
# ----------------------------------------
#  based on tcshell scritps from prof Roel Neggers
#  scripted and checked by Jan Chylik, August 2016 
#  inquiries and error reports: jchylik@uni-koeln.de
#  
# last update: 21 September 2016, 11 January 2017

# === setting  ================================================
# if the script is run from dales_prepare.sh, the setting does not require any changes

# ---- directory setting ---------------------------------
#  directroy in which to run 
    export rundir=./
    # # export rundir=/work/${USERN}/dales/dales4/experiments/Juelich/juelich20130605_imicro2_16x16
# busy file 
    export busyfile="runDALES_busy.txt"
# runlog
    export RUNLOG="runlog.txt"
# setting 
    export N_MPI=1
    export DALES="dales4"
# === script ==============================================
 
#------------- important directories ----------------------


   
# -------- run  itself ------------------------------------

if [ ! -e $busyfile ]; then

  #--- create file indicating that job is ongoing ---
  echo "runDALES job is busy ... don't interfere!" > $busyfile

  cd ${rundir}

  #echo "test" > test.txt

  # running dales
  echo "Starting dales run, you can view the run log  ${rundir}$RUNLOG "
  echo "You can view it in other terminal."
  echo "Do not close this terminal until the job is finished."
  #./dales_roel4_20140423
  # mpirun -np $N_MPI ./${DALES_VERSION}_${USERN}_${DAYC}_${MACHINE_TYPE}>$RUNLOG
  mpirun -np $N_MPI ${DALES}>$RUNLOG
  # note: it keeps the name of the original executable so we can have more different dales executables possibly used alongside each other if needed

  #--- cleanup --- 
  rm $busyfile
  
  echo " "
  echo "Dales run has ended."
  exit 0

else

  echo " "
  echo "A runDALES job is already in progress ... aborting"
  echo "Exitting with exit status 3"
  exit 3

fi


         
# ===========================================================


