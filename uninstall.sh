#!/bin/bash

# Removes files/folders created during a CMake build
# Run with:
# ./uninstall


for i in `cat build/install_manifest.txt`
do
	echo "Removing installed file $i"
	rm -f $i
done

echo "Removing directory build"

rm -rf build
