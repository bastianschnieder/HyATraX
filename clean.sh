#!/bin/bash

for dir in */; do
    if [ -e "$dir/control" ]; then
        (cd "$dir"
        mv BACK/* .
	rm -r BACK
        rm err* out* define_level2.in define.in)
    fi
done
