echo "Expected output:"
grep "Elapsed" expected_output/output
echo "Current output:"
grep "Elapsed" output
echo " "

len1=`wc expected_output/path.info | awk '{print $1}'`
len2=`wc path.info | awk '{print $1}'`

if [ $len1 != $len2 ]; then
  echo "path.info files are of different lengths. Test failed."
  exit 1
fi

python path_test_parser.py
exit $?
