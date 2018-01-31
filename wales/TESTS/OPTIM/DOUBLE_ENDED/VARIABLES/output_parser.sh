echo "Expected output:"
grep "Elapsed" expected_output/output
echo "Current output:"
grep "Elapsed" output
echo " "

natoms1=`grep "getparams> Number of atoms" expected_output/output | awk '{print $9}'`
natoms2=`grep "getparams> Number of atoms" output | awk '{print $9}'`

if [ $natoms1 != $natoms2 ]; then
  echo "Different numbers of optimisation coordinates in the systems. Test failed."
  exit 1
fi

len1=`wc expected_output/path.info | awk '{print $1}'`
len2=`wc path.info | awk '{print $1}'`

if [ $len1 != $len2 ]; then
  echo "path.info files are of different lengths. Test failed."
  exit 1
fi

ntriples1=`grep "TS in the path" expected_output/output | awk '{print $8}'`
ntriples2=`grep "TS in the path" output | awk '{print $8}'`

if [ $ntriples1 != $ntriples2 ]; then
  echo "Different discrete path lengths. Test failed."
  exit 1
fi

python connect_test_parser.py $ntriples1 $natoms1
exit $?
