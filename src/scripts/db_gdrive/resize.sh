#!/bin/bash

for file in imgs/*.png; do
    convert $file -resize 100x100 $file;
done
