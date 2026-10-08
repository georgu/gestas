#!/bin/sh
#
# compares the source files with their originals
#
#-------------------------------------------------------

compdir="$HOME/shyfem/src/contrib/ecological/weutro"
utildir="$HOME/shyfem/src/utils/shyutilmpi"

files=$( ls *.f90 )

for file in $files
do
  echo $file
  dir=$compdir
  [ $file = "shympi_fake.f90" ] && continue
  [ $file = "writevars.f90" ] && dir=$utildir
  cmp $dir/$file $file
  if [ $? -ne 0 ]; then
    tkdiff $dir/$file $file
  fi
done

#-------------------------------------------------------

