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

if [ ! -f $file ]; then
  echo "*** no such file: $file"
  exit 1
fi

rm -f box_*.tmp
splitcol $file

for col in col.*
do
  number=$( echo $col | sed -e 's/col.//' )
  newnumber=$number
  [ $number -lt 10 ] && newnumber="0$newnumber"
  echo "$col $number -> $newnumber"
  mv $col box_$newnumber.tmp
done

title=$( head -1 $file )
gp -c -t "file $file - $title" -tx date -ty value -l ".tmp" box_*.tmp

gv out.ps

#------------------------------------------------

