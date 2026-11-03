#!/bin/sh
#
# converts simulation files from the old to the new convention
#
#--------------------------------------------------------------

if [ $# -ne 1 ]; then
  echo "Usage: convertsim.sh simname"
  exit 1
fi

simname=$1

CheckFile()
{
  local what=$1
  local file="boxes_$what.txt"

  if [ ! -f $file ]; then
    echo "No such file: $file ...aborting"
    exit 1
  fi
}

RenameFile()
{
  local what=$1
  local file="boxes_$what.txt"

  if [ -f $file ]; then
    newfile="$simname.boxes.$what.txt"
    echo "  $file  ->  $newfile"
    mv $file $newfile
  fi
}

#--------------------------------------------------------------

CheckFile 2d
CheckFile geom
CheckFile init
CheckFile meteo
CheckFile stats

#--------------------------------------------------------------

RenameFile 2d
RenameFile geom
RenameFile init
RenameFile meteo
RenameFile stats

RenameFile 3d
RenameFile vertical
RenameFile matrix

#--------------------------------------------------------------

