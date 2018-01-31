#!/bin/sh
#
# Written by Ilyas Yildirim, Dept. of Chem and Biochem, Florida Atlantic University (March 2, 2017)
#
# This script manipulates the .ps file produced by disconnectionDPS. The scales are optimized for my own purposes. Change it the way you want.
#

scale=4

cat tree.ps | \
awk -v sc=$scale '{ \
  if((/ mt /) && (/ ls$/)){\
    $1=$1*sc; \
    $4=$4*sc; \
    fmt++; \
    if(fmt == 1){\
      refalign= $1; \
    } \
  }else if((/ mt /) &&(/show/)){\
    $1=$1*sc; \
    if(f==1){ \
      $1=refalign-70; \
    } \
  } else if(/scalefont/) {\
    f++; \
    $3=$3*sc/2; \
  } else if(/BoundingBox/) {\
    $4=$4*sc; \
  } \
  print $0; \
}' > tree.mod.ps
