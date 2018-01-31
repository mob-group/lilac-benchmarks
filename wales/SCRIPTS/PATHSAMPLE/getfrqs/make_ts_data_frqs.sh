awk '{printf "%20s\n", $1}' ts.data.removed > energy.tmp
awk '{printf "%20s\n", $2}' frqs.sorted.ts > frqs.tmp
paste  energy.tmp frqs.tmp> ts.data.frqs.tmp
awk '{printf "%4s %6s %6s %20s %20s %20s\n", $3, $4, $5, $6, $7, $8}' ts.data.removed > tmp
paste   ts.data.frqs.tmp tmp > ts.data.frqs
rm *tmp
