#!/bin/bash

# Set default workin dir, if no second argument is given
if [ -z "$2" ]; then
  WDIR=/Users/rei/lab/AGN_relativistic
else
  WDIR=$2
fi

# Link definitions-$1.h if no link exists yet
if [ ! -e definitions.h ]; then
  ln -s definitions-$1.h definitions.h
  echo "linked definitions.h to definitions-$1.h"

# If link exists and link and chosen definitions-$1.h are different,
# remove link first, and the relink.
elif [ "$(diff definitions-$1.h definitions.h)" ]; then
  rm definitions.h
  ln -s definitions-$1.h definitions.h
  echo "relinked definitions.h to definitions-$1.h"

# If link and chosen definitions-$1.h are the same, do not do anything
fi

cp pluto-$1.ini $WDIR/$1/pluto.ini
cp pluto $WDIR/$1/
cp sysconf.out $WDIR/$1/
cp cooltable.dat $WDIR/$1/

# cp input-rho_256_12.flt $WDIR/$1/
# cp grid_in.out $WDIR/$1/

