for f in $@
do
    f1=$f\.dat
    echo gawk -f convert.gawk $f \> $f1
    gawk -f convert.gawk $f > $f1
done 
