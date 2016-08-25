#!/bin/bash
for file in $@
do
    out=$file\.tmp
    gawk -f rewrite.gawk $file > $out
    echo $out
    mv $out $file
done
