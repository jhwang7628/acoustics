#/bin/bash

rename 's/\_/\./g' head_*
rename 's/\./\_/' head.*
echo ">> striping original headers"
ls head_* | awk '{print "cat " $1 " | tail -n +2 > " $1 "--" }' | sh
rm head_*.?????
rename 's/\-\-//' head_*

echo ">> pasting positions"
ls head_* | awk '{print "paste -d, position.txt " $1 " > " $1 "--" }' | sh
rm head_*.?????
rename 's/\-\-//' head_*
rename 's/\./\.csv\./' head_*

