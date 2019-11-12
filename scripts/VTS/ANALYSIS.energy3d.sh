#!/bin/bash
# submit evndisp energy3d
#

# qsub parameters
h_cpu=11:59:00; h_vmem=12000M; tmpdir_size=15G
#h_cpu=23:59:00; h_vmem=20000M; tmpdir_size=25G

if [ $# -lt 5 ]; then
# begin help message
echo "
--------------------------------------------------------------------------------
Energy3D adds energy to 3D-modeled files

ANALYSIS.energy3d.sh <REAL/SIM>> <Runlist> <Noicelist> <Filedir> <Refdir> <Outdir> <CAREGRISU> <TelVerlist> <SumWintlist>

required parameters:

	<REAL/SIM>	Runlist with info if the run is real or simulated data, REAL/SIM

	<Runlist>	Runlist with files that have been 3D-modeled, .txt file with runnumbers 

	<Noicelist>	List with the noiceleves for the runlist (use pedvarckeck to get this)

	<Filedir>	List with the directorys where the 3D-modeled files are, ''model3d'' for REAL files, zenithangle for SIM (ex. ''35'')

	<Refdir>	The directory where the 3D-model reference files are

	<Outdir>	Your wanted output directory

	<CAREGRISU>	Runlist with simulation choises, CARE/GRISU

	<TelVerlist> 	Runlist with the versions of telescope setup, 4/5/6

	<SumWintlist> 	Runlist with winter or summer atmosphere, 21/22
	
Example: ./ANALYSIS.energy3d.sh listdir/RealSimList.txt listdir/RunList.txt listdir/NoiceList.txt listdir/DirectoryList.txt $VERITAS_USER_DATA_DIR/DATASets3DFull/Merged/ $VERITAS_USER_DATA_DIR/analysis/energy3d listdir/caregrisulist.txt listdir/telverlist.txt listdir/sumwintlist.txt 
--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Run init script
bash "$( cd "$( dirname "$0" )" && pwd )/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# EventDisplay version
EDVERSION=`$EVNDISPSYS/bin/evndisp --version | tr -d .`

# Parse command line arguments
REALORSIMLIST=$1	#list
RUNLIST=$2		#list
NOISELIST=$3		#list
DATADIRLIST=$4		#list
REFDIR=$5		#not list
OUTDIR=$6		#not list
CAREGRISULIST=$7	#list
TELVERLIST=$8		#list
SUMWINTLIST=$9		#list

# directory for run scripts
DATE=`date +"%y%m%d"`
LOGDIR="$VERITAS_USER_LOG_DIR/$DATE/EVNDISP.ANAMCVBF"
mkdir -p $LOGDIR

# output directory for evndisp products (will be manipulated more later in the script)
if [[ ! -z "$OUTDIR" ]]; then
    ODIR="$OUTDIR"
fi
# output dir
OPDIR=$OUTDIR
mkdir -p $OPDIR
chmod -R g+w $OPDIR
echo -e "energy3d output files will be written to:\n $OPDIR"
  
# Job submission script
SUBSCRIPT="$EVNDISPSYS/scripts/VTS/helper_scripts/ANALYSIS.energy3d_sub"
  
# reading from the 7 lists
RSFILS=`cat $REALORSIMLIST`
NREALSIMS=`cat $REALORSIMLIST | wc -l `

RFILES=`cat $RUNLIST`
NRUNS=`cat $RUNLIST | wc -l `

NFILES=`cat $NOISELIST`
NNOISES=`cat $NOISELIST | wc -l `

DATVFILES=`cat $DATADIRLIST`
NDIRS=`cat $DATADIRLIST | wc -l  `

TELFILES=`cat $TELVERLIST`
NTELVERS=`cat $TELVERLIST | wc -l  `

SWFILES=`cat $SUMWINTLIST`
NSUMWINS=`cat $SUMWINTLIST | wc -l  `

CGFILES=`cat $CAREGRISULIST`
NCAREGRISUS=`cat $CAREGRISULIST | wc -l  `



# loop over each run
h=1;
for RSFIL in $RSFILS; do
  i=1;
  for RFIL in $RFILES ; do
    j=1;
    for NFIL in $NFILES ; do 
      k=1;
      for DFIL in $DATVFILES ; do 
        l=1;
        for TVFIL in $TELFILES ; do 
          m=1;
          for SWFIL in $SWFILES ; do 
            n=1;
              for CGFIL in $CGFILES ; do 
                 if [ $h -eq $i -a  $i -eq $j -a  $j -eq $k  -a $k -eq $l -a $l -eq $m -a $m -eq $n ] ; then 
		   if [ $RSFIL == "REAL" ] ; then	
                      echo -e "Run being submitted: $RFIL with noice $NFIL"
		   else		
                     echo -e "Run being submitted: $RFIL with noice $NFIL at $DFIL dergees zenithangle"
   		   fi
                   # make run script
                   FSCRIPT="$LOGDIR/EVN.energy3d_$RFIL"
    
                   echo -e "ALL VARIABLES SENDING $RSFIL $RFIL $NFIL $DFIL $REFDIR $OUTDIR $TVFIL $SWFIL $CGFIL"
  
                   sed -e "s|REALSIM|$RSFIL|"\
		       -e "s|RNUM|$RFIL|" \
                       -e "s|NOICELEVELS|$NFIL|" \
                       -e "s|INDATADIR|$DFIL|" \
                       -e "s|RDIR|$REFDIR|" \
                       -e "s|OUTPUTDIR|$OUTDIR|" \
                       -e "s|TELV|$TVFIL|" \
                       -e "s|SUMW|$SWFIL|" \
                       -e "s|CAREGRISU|$CGFIL|" $SUBSCRIPT.sh > $FSCRIPT.sh  
	
                   chmod u+x $FSCRIPT.sh
                   echo $FSCRIPT.sh

                   # run locally or on cluster
                   SUBC=`$EVNDISPSYS/scripts/VTS/helper_scripts/UTILITY.readSubmissionCommand.sh`
                   SUBC=`eval "echo \"$SUBC\""`
                   if [[ $SUBC == *"ERROR"* ]]; then
                      echo $SUBC
                      exit
                   fi  
       
                   if [[ $SUBC == *qsub* ]]; then
                     JOBID=`$SUBC $FSCRIPT.sh`      
                     echo "RUN $RFIL: JOBID $JOBID"
                   elif [[ $SUBC == *parallel* ]]; then
                     echo "$FSCRIPT.sh &> $FSCRIPT.log" >> $LOGDIR/runscripts.dat
                   fi
                 fi
               ((n++))
             done
             ((m++))
           done
           ((l++))
         done
         ((k++))
       done
       ((j++))
     done
     ((i++))
   done
   ((h++))
done
                               
exit
