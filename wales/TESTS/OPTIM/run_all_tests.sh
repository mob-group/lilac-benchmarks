exec=${1}

if [ -z "${exec}" ]; then
  echo "No executable specified."
  exit 1
fi

bash clean_directories.sh
touch tests_output
echo "New exec: ${exec}" >> tests_output


for test in GEOPT PATH DOUBLE_ENDED FRQS; do

#    for sys in AMBER12 AMBER9 AMBER9_genrigid BLJ60 BOWMAN LJ38 SD TIP4P_genrigid TIP4P_rbaa TTM3 VARIABLES; do
    for sys in AMBER12 BLJ60 LJ38 SD TIP4P_genrigid TIP4P_rbaa TTM3; do  # Subset of tests to run if you're in a hurry!

	echo "********************Running test ${test} for system ${sys}********************"
	echo "********************Running test ${test} for system ${sys}********************" >> tests_output

	cd ${test}/${sys}

	if [ "${test}" == "FRQS" ]; then
	    for test2 in GEOPT PATH DOUBLE_ENDED; do
		cd ${test2}
		${exec} > output
		cd ..
	    done
        else
	    ${exec} > output
	fi

	bash output_parser.sh >> ../../tests_output 2>&1
	cd ../..

     done
done
  
