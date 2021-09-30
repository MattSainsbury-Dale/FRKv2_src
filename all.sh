#!/bin/bash

set -e

echo ""
echo "######## SETUP ############"
echo ""

## Should we use "quick", very-low-dimensional representations of the models?
echo "Do you wish to use very-low-dimensional representations of the models to quickly establish that the code is working (note that the generated results and plots will not exactly match those in the manuscript if you reply 'y', and you will need at least 32GB of RAM and/or swap space if you reply 'n')? (y/n) "
read quick_str

if ! [[ $quick_str == "y" ||  $quick_str == "Y" || $quick_str == "n" || $quick_str == "N" ]]; then
    echo "Please re-run and type type y or n"
    exit 1
fi


if [ $quick_str == "y" ] || [ $quick_str == "Y" ]; then
    quick=1
else
    quick=0
fi
## echo quick = $quick

Rscript setup.R

echo ""
echo "######## STARTING POISSON EXAMPLE OF SECTION 3.1 ############"
echo ""
Rscript scripts/Poisson_sim.R $quick
echo ""

echo "######## STARTING NEGATIVE-BINOMIAL EXAMPLE OF SECTION 3.2 #############"
echo ""
Rscript scripts/Negbinom_sim.R $quick
echo ""

echo "######## STARTING HEATON COMPARISON OF SECTION 3.3 #############"
echo ""
Rscript scripts/Heaton.R $quick
echo ""

echo "######## STARTING MODIS COMPARISON OF SECTION 4.1 #############"
echo ""
Rscript scripts/MODIS.R $quick
echo ""

echo "######## STARTING AMERICIUM COMPARISON OF SECTION 4.2 #############"
echo ""
Rscript scripts/Am.R $quick
echo ""

echo "######## STARTING SYDNEY EXAMPLE OF SECTION 4.3 #############"
echo ""
Rscript scripts/Sydney.R $quick
echo ""

echo "######## STARTING CHICAGO EXAMPLE OF SECTION 4.4 #############"
echo ""
Rscript scripts/Chicago.R $quick
echo ""
