#!/bin/bash
for dir in */; do
    if [ -e "$dir/control" ]; then
        (cd "$dir"
        mv BACK/* .
        rm -r BACK err* out*
        cp /home/USERNAME/ts_opt.sh .
        cp /home/USERNAME/define_level2.in .
        define < define_level2.in
        kdg end
        echo "\$statpt" >> control
        echo "   itrvec      1" >> control
	echo "   tradius     3.000000E-02" >> control
        echo "\$last step    define" >> control
	echo "\$arh" >> control
        echo "\$end" >> control
        sbatch -J "$dir" ts_opt.sh
        echo "done for $dir")
    fi
done
