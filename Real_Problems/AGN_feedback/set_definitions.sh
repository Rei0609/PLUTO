#!/bin/bash

WDIR=/Users/rei/lab/AGN_feedback

if [ ! -e definitions.h ]; then
  ln -s definitions-$1.h definitions.h
  echo "linked definitions.h to definitions-$1.h"

elif [ "$(diff definitions-$1.h definitions.h)" ]; then
  rm definitions.h
  ln -s definitions-$1.h definitions.h
  echo "relinked definitions.h to definitions-$1.h"

fi

cp pluto-$1.ini $WDIR/$1/pluto.ini
cp pluto $WDIR/$1/
cp sysconf.out $WDIR/$1/
# cp input-rho_256_12.flt $WDIR/$1/
# cp grid_in.out $WDIR/$1/

