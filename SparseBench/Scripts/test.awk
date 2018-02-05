BEGIN {validate=0; }
/Iteration/ {if ($5>1.e-12) {iter=$2; validate=$5}}
/Matrix/ || /Preconditioner/ || /Iterative/ || /Vector/ || /Overall/ {print}
/Mflop/ {print}
/Dimension/ {print}
/Integer:/ {print}
/Real:/ {print}
END {print "Final error: ",iter,validate; print ""}
