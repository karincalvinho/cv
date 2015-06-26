#!/bin/sh
echo "Filename,Tafel Slope (mV/dec),Potential at 1 mA/cm2 vs NHE (V),Potential at 10mA/cm2 (V),Exchange Current Density (mA/cm2),Intercept,R-value"
find ./data/02-2015 -name "CV*2.txt" | while read x
do
  path=$(dirname "${x}")
  filename=$(basename "${x}")
  python cvs.py "$x"
done
find ./data/03-2015 -name "CV*2.txt" | while read x
do
  path=$(dirname "${x}")
  filename=$(basename "${x}")
  python cvs.py "$x"
done
find ./data/04-2015 -name "CV*2.txt" | while read x
do
  path=$(dirname "${x}")
  filename=$(basename "${x}")
  python cvs.py "$x"
done
find ./data/05-2015 -name "CV*2.txt" | while read x
do
  path=$(dirname "${x}")
  filename=$(basename "${x}")
  python cvs.py "$x"
done
find ./data/06-2015 -name "CV*2.txt" | while read x
do
  path=$(dirname "${x}")
  filename=$(basename "${x}")
  python cvs.py -c ch_145.yaml "$x"
done
