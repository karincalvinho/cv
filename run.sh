#!/bin/sh
echo "Filename,Tafel Slope (mV/dec),Potential at 0.1 mA/cm2 vs NHE (V),Potential at 10mA/cm2 (V),Exchange Current Density (mA/cm2),Intercept,R-value"
find ./data/March2015/033015 -name "CV*.txt" | while read x
do
  path=$(dirname "${x}")
  filename=$(basename "${x}")
  python cvs.py "$x"
done
