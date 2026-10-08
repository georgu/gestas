#!/bin/sh
#
# plots results of boxes
#
#------------------------------------------------

if [ $# -lt 1 ]; then
  echo "Usage: plotboxes.sh result-file"
  exit 1
fi

file=$1

rm -f box_*.tmp
splitcol $file

for col in col.*
do
  number=$( echo $col | sed -e 's/col.//' )
  #echo "$col $number"
  newnumber=$number
  [ $number -lt 10 ] && newnumber="0$newnumber"
  echo "$col $number -> $newnumber"
  mv $col box_$newnumber.tmp
done

gp -c -l ".tmp" box_*.tmp

gv out.ps

#------------------------------------------------

