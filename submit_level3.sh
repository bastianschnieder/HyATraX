#!/bin/bash

for dir in */; do
    if [ -e "$dir/control" ]; then
        (cd "$dir"
        mv BACK/* .
        rm -r BACK 
	cp /home/USERNAME/submit.sh .
        sbatch -J "$dir" submit.sh)
    fi
done
