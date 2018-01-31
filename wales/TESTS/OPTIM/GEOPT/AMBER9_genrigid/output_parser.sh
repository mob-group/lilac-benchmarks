echo "Expected output:"
grep "Elapsed" expected_output/output
echo "Current output:"
grep "Elapsed" output
echo " "

python geopt_test_parser.py
exit $?
