#!/bin/sh

maxmin=15673927

cp min.data min.data.$maxmin

for sigma in 0.02 0.03 0.04 0.05 0.06 0.07 0.08 ; do
   rm structures.$sigma
   echo 1 32 1000 1000 0.04 0.375  500 219 na    na    all 0.01 $sigma 10.0 > Cv.data
   cat Cv.data.template >> Cv.data

   for nmin in 10000 15000 25000 50000 100000 1000000 2000000 3000000 4000000 5000000 6000000 7000000 8000000 10000000 12000000 14000000 $maxmin ; do

       rm merge
#      for mergei in 1 ; do
#          for mergeq in 6 8 10 12 14 16 18 20 ; do

               head -$nmin min.data.$maxmin > min.data
#              echo $mergei > merge
#              echo $mergeq >> merge
               ~/svn/UTILS/GMIN/calc.Cv > Cv.output.$nmin.$sigma
               cp minima.pdf minima.pdf.$nmin.$sigma

               echo $nmin `grep "Total number of structures" Cv.output.$nmin.$sigma \
                    | sed -e 's/Total number of structures is//'` >> structures.$sigma
                           
#          done
#      done
   done
done
