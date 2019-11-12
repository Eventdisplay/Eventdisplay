# script to upload packages for releases to UCLA

UCLA=`grep "\* VTSRAWDATA" $VERITAS_EVNDISP_AUX_DIR/ParameterFiles/EVNDISP.global.runparameter | awk '{print $3}'`

echo $UCLA

# upload version might be different from download version
UP="570"
VERSION="EVNDISP-${UP}-auxv01"

bbftp -u bbftp -V -S -m -p 12 -e "put EVNDISP-${UP}.tar.gz /veritas/upload/EVNDISP/EVNDISPDATA.v${UP}/EVNDISP-${UP}.tar.gz" $UCLA

# for N in lookuptables runfiles dispBDTs Model3D effectiveareas radialacceptances Frogs
# for N in lookuptables runfiles radialacceptances calibration Frogs dispBDTs Model3D
# for N in effectiveareas.V5-ATM21 effectiveareas.V5-ATM22 effectiveareas.V4-ATM21 effectiveareas.V4-ATM22
for N in dispBDTs
do
   FILE=${VERSION}.VTS.aux.${N}.tar.gz
   echo $FILE
   ls -lh $FILE
   md5sum $FILE
   bbftp -u bbftp -V -S -m -p 12 -e "put $FILE /veritas/upload/EVNDISP/EVNDISPDATA.v$UP/$FILE" $UCLA
done
