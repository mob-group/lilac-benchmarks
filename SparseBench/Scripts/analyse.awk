/Dimension/ {d=$2; mem=0;}
/Integer/ || /Real/ {mem=mem+$4;}
/Matrix/ {m=1}         /Mflop/ && m==1 {m=$3}
/Preconditioner/ {p=1; i=0} /Mflop/ && p==1 {p=$3}
/Iterative/ {i=1} /Mflop/ && i==1 {i=$3}
/Vector/ {v=1}         /Mflop/ && v==1 {v=$3}
/Overall/ {o=1}
/Mflop/ && o==1 {
    if (m<.001 || m>1.e6) m=0;
    if (p<.001 || p>1.e6) p=0;
    if (v<.001 || v>1.e6) v=0;
    printf("Size: %d Mb= %7.4f ; MFlop rate: %7.3f m= %7.3f p= %7.3f v= %9.3f\n",d,mem,$3,m,p,v);
    o=0;}
