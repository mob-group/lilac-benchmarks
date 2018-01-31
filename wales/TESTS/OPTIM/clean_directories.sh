
for test in GEOPT PATH DOUBLE_ENDED FRQS; do

#    for sys in AMBER12 AMBER9 AMBER9_genrigid BLJ60 BOWMAN LJ38 SD TIP4P_genrigid TIP4P_rbaa TTM3 VARIABLES; do
    for sys in AMBER12 BLJ60 LJ38 SD TIP4P_genrigid TIP4P_rbaa TTM3 VARIABLES; do  # Quick list!

	cd ${test}/${sys}

	if [ "${test}" == "FRQS" ]; then
	    for test2 in GEOPT PATH DOUBLE_ENDED; do
		cd ${test2}
		rm min.data.info output points* energies odata.new chirality_readable initial_c*_states min.out EofS* path.* fort.* neb* *ofI rbneb* rbpath* rbsites* > /dev/null 2>&1
		cd ..
	    done
        else
	    rm min.data.info output points* energies odata.new chirality_readable initial_c*_states min.out EofS* path.* fort.* neb* *ofI rbneb* rbpath* rbsites* > /dev/null 2>&1
	fi

	cd ../..

     done
done

  
