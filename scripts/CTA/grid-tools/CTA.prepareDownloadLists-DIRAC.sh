# script to prepare donwload lists

if [ ! -n "$1" ]  && [ ! -n "$2" ]
then
   echo "./CTA.prepareDownloadLists-DIRAC.sh <site> <set>"
   echo
   echo ".e.g for Paranal_proton_South_20deg_HB9 do: "
   echo "./CTA.prepareDownloadLists-DIRAC.sh Paranal_ _20deg_HB9"
   exit
fi

# DIRAC variables must be set
if [ -z "$DIRAC" ]
then
   echo "Error: set dirac environement"
   exit
fi
   
SITE="$1"
SET="$2"

# Paranal_proton_South_20deg_HB9

for P in gamma gamma-diffuse proton electron
do
   for A in North South
   do
       DSET="${SITE}${P}_${A}${SET}"
       echo $DSET
       cta-prod3-dump-dataset $DSET
   done
done
