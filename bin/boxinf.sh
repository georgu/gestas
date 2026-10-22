#!/bin/sh

dir=$HOME/work/boxes/estas

posts="geom stats init 2d meteo"

actdir=$( pwd )
cd $dir
make boxinf
[ $? -ne 0 ] && echo "errors compiling..." && exit 1
cd $actdir

for post in $posts
do
  file=boxes_$post.txt
  [ ! -f $file ] && continue
  echo "======================================="
  echo $file
  echo "======================================="
  $dir/boxinf < $file
done

