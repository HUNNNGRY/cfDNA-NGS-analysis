#!/bin/bash 
## perform multi join func. in linux (1,2  +  1,3  -->  1,2,3)
out=$1
shift 1
cat $1 | awk '{print $1}' > $out
for f in $*; do join -1 1 $out -2 1 $f > tmp1; mv tmp1  $out; done

