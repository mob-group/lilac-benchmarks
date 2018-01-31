#!/bin/sh

cp pathdata pathdata.save
echo "EXTRACTMIN 1" > pathdata
cat pathdata.save >> pathdata
/sharedscratch/wales/PATHSAMPLE.pgi/PATHSAMPLE
nvar=`wc extractedmin | awk '{print $1}'`
echo number of variables is $nvar

echo "EXTRACTMIN -1" > pathdata
cat pathdata.save >> pathdata
/sharedscratch/wales/PATHSAMPLE.pgi/PATHSAMPLE

cp pathdata.save pathdata
rm AUC.test
rm AUC.pasted
count=-1

for start in 1 1000 1500 2000 2500 ; do

   string=`head -1 odata.multi.prob`
   echo $string $start > odata
   sed -e "1d" odata.multi.prob >> odata
   head -$nvar extractedmin >> odata
   /sharedscratch/wales/OPTIM.pgi/OPTIM > output.multi.prob.$start

   grep AUC output.multi.prob.$start | sort -u | awk '{print $4}' > AUC.$start
   cat AUC.$start >> AUC.test
   wc min.data | awk '{print $1}'> poo
   mv poo min.AUC.$start
   cat  AUC.test >> min.AUC.$start
   if [ $count == -1 ]; then
      cp AUC.$start AUC.pasted
   else
      paste AUC.pasted AUC.$start > temp
      cp temp AUC.pasted
   fi
   ((count++))

done

echo $count `wc AUC.1` > AUCdata
cat AUC.test >> AUCdata
../../../../AUCmean < AUCdata > AUC.mean
paste AUC.pasted AUC.mean > temp
cp temp AUC.pasted
gmin=`cat -n min.data | sort -k 1n | head -1 | awk '{print $1}'`
sed -n ${gmin}p AUC.pasted > AUC.gmin

