#!/bin/bash

# Takes 'test' and 'control' directories as arguments. 
DIR1=$1
DIR2=$2

# The directories should contain subdirectories with fits files which match names
# between the 'test' and 'control' directories.
# e.g `test/batch1/name1.fits` matches to `control/batch1/name1.fits`

DIFF_FOUND=false
# Use `echo $?` to retrieve whether the script found any differences 
# (0 = no differences, 1 = at least one file contained differences)


# Works recursively
for FILE in $(find $DIR1 -name \*.fits); do
    fitsdiff $FILE ${FILE/$DIR1/$DIR2} > $FILE.txt
    
    if ! grep -q "No differences found" $FILE.txt; then
        cat $FILE.txt
        DIFF_FOUND=true
    fi
done

if [ "$DIFF_FOUND" = true ] ; then
    exit 1
fi