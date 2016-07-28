#!/bin/bash 

# remove all files and directories from the current location.
for p in *; do 
    echo "Removing $p.."
    rm -r $p
done
