#!/bin/bash

WDIR=/home/rei/research/lab/AGN_wind/sim

if [ ! -e definitions.h ]; then
  ln -s definitions-$1.h definitions.h
  echo "linked definitions.h to definitions-$1.h"

elif [ "$(diff definitions-$1.h definitions.h)" ]; then
  rm definitions.h
  ln -s definitions-$1.h definitions.h
  echo "relinked definitions.h to definitions-$1.h"

fi

cp pluto-$1.ini $WDIR/$1/pluto.ini

