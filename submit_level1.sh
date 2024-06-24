#!/bin/bash

for file in coord_hidx_*.xyz;
do
	num="${file#coord_}"
	mkdir -p "$num"
	mv "$file" "$num/"
done
for dir in */; do
        (cd "$dir"
        cp /home/USERNAME/submit.sh .
        x2t *.xyz > coord
	cp /home/USERNAME/define.in .
        define < define.in)
done

for dir in */; do
    if [ -e "$dir/control" ]; then
        (cd "$dir"
        kdg end
        kdg scfiterlimit
        echo "\$scfiterlimit 800" >> control
        kdg maxcor
        echo "\$maxcor    3000 MiB  per_core" >> control
        kdg dft
        echo "\$dft" >> control
        echo "   functional tpss" >> control
        echo "   gridsize   m3" >> control
        echo "   weight derivatives" >> control
        echo "\$end" >> control
        rm define.in
        sbatch -J "$dir" submit.sh
        echo "done for $dir")
    fi
done


