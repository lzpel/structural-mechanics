#!/bin/sh
set -eux
gfortran $1
chmod 777 ./a.exe
cat $2 | ./a.exe
echo "set size ratio -1 ; unset border ; unset tics ; unset key ; a=10. ; plot 'deformation.dat' using (\$1+a*\$3):(\$2+a*\$4) with lines" | gnuplot -p -
