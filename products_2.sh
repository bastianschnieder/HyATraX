#!/bin/bash
for dir in */; do
    if [ -e "$dir/control" ]; then
        (cd "$dir"
        mv BACK/* .
        rm -r BACK err* out*
	cp /home/USERNAME/submit.sh .
        sbatch -J "$dir" submit.sh
        echo "done for $dir")
    fi
done
