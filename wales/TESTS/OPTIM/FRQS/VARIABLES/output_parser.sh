natoms1=`grep "getparams> Number of atoms" DOUBLE_ENDED/expected_output/output | awk '{print $9}'`
natoms2=`grep "getparams> Number of atoms" DOUBLE_ENDED/output | awk '{print $9}'`

if [ $natoms1 != $natoms2 ]; then
  echo "Different numbers of optimisation coordinates in the systems. Test failed."
  exit 1
fi

len1=`wc DOUBLE_ENDED/expected_output/path.info | awk '{print $1}'`
len2=`wc DOUBLE_ENDED/path.info | awk '{print $1}'`

if [ $len1 != $len2 ]; then
  echo "path.info files are of different lengths. Test failed."
  exit 1
fi

nTS1=`grep "TS in the path" DOUBLE_ENDED/expected_output/output | awk '{print $8}'`
nTS2=`grep "TS in the path" DOUBLE_ENDED/output | awk '{print $8}'`

if [ $nTS1 != $nTS2 ]; then
  echo "Different discrete path lengths. Test failed."
  exit 1
fi

python frqs_tests_parser.py $nTS1 $natoms1
exit $?
