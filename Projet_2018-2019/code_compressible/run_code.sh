#!/bin/bash

rm -f resultat* *.pdf

echo 'run the Euler Equation with HLL solver'
cd Schema_HLL
make clean && make && ./run
cp resultat_HLL.txt ../
cd ..

echo 'run the Euler Equation with HLLC solver'
cd Schema_HLLC
make clean && make && ./run
cp resultat_HLLC.txt ../
cd ..

echo 'plot les resultats'
gnuplot  affichage_HLLC_HLL.gnu
