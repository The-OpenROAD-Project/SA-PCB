#!/bin/bash
i=5000000
j=20
t=40000

path="./benchmarks/mcnc/"
cname="apte"

for k in {1..5}
do
  plfname="${cname}_thread_${k}.pl"
   ./sa -i $i/$k -j $j -t $t -k $k -s -w 0 -f $plfname
   echo "======"
   python evaluate "$path$apte.nodes $plfname $path$apte.nets"
   python plot "$path$apte.nodes $plfname $path$apte.nets"
   echo "======""
done
