awk '{printf "%20s\n", $1}' min.data.removed > energy.tmp
awk '{printf "%20s\n", $2}' frqs.sorted.min > frqs.tmp
paste  energy.tmp frqs.tmp> min.data.frqs.tmp
awk '{printf "%4s %20s %20s %20s\n", $3, $4, $5, $6}' min.data.removed > tmp
paste   min.data.frqs.tmp tmp > min.data.frqs
rm *tmp
