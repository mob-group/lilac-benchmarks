/Dimension/ {d=$2; mem=0;}
/Integer/ || /Real/ {mem=mem+$4;}
/Matrix/ {m=$4}
/Preconditioner/ {p=$4}
/Vector/ {v=$4}
/Overall/ {
    if (m<.001 || m>1.e4) m=0;
    if (p<.001 || p>1.e4) p=0;
    if (v<.001 || v>1.e4) v=0;
    printf("Size: %d Mb= %7.4f MFlop rate: %7.3f m=%7.3f p=%7.3f v=%7.3f\n",d,mem,$4,m,p,v);}
