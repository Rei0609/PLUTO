#!/bin/bash

WDIR=/Users/ayw/lab/Radiative_Shocks_in_Clouds/sim

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
cp cooltable.dat $WDIR/$1/
cp cooltable_frac.dat $WDIR/$1/

