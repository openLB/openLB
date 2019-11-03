#!/bin/bash

# Script to run the example with three different values for maxU and collect the output.

APP=powerLaw2d
NP=2

if [ ! -f $APP ]; then make; fi

for N in `seq 1 3`; do

		  case $N in
					 1) MAXU=1;   DIR="out_u1e+0" ;;
					 2) MAXU=0.1; DIR="out_u1e-1" ;;
					 3) MAXU=10;  DIR="out_u1e+1" ;;
		  esac

		  if [ ! -d $DIR ]; then mkdir $DIR; fi
		  
		  sed "s/MAXU/$MAXU/g" input-model.xml > $DIR/input.xml
		  cp $APP $DIR/

		  cd $DIR

		  mpirun -np $NP $APP 2>&1 | tee log.log \
					 && cp tmp/gnuplotData/centerVelocity.png ../u_$DIR".png"

		  cd ..
done
