@echo off

echo "######## SETUP ############"

:choice
set /P Input=Do you wish to use very-low-dimensional representations of the models to quickly establish that the code is working (note that the generated results and plots will not exactly match those in the manuscript if you reply 'Y', and you will need at least 32GB of RAM and\or swap space if you reply 'N')[Y/N]?
if /I "%Input%" EQU "Y" goto :quickchoice
if /I "%Input%" EQU "N" goto :nonquickchoice
goto :choice


:quickchoice
set /A quick = 1
echo "Run with very-low-dimensional representations of the models, quick = %quick%"
goto :callRscript

:nonquickchoice
set /A quick = 0
echo "Not run with very-low-dimensional representations of the models, quick =%quick%"
goto :callRscript

:callRscript
Rscript scripts\setup.R


echo "######## STARTING SYDNEY EXAMPLE OF SECTION 3.1 #############"
Rscript scripts\Sydney.R %quick%

echo


echo "######## STARTING CHICAGO EXAMPLE OF SECTION 3.2 #############"
Rscript scripts\Chicago.R %quick%
echo


echo "######## STARTING HEATON COMPARISON OF SECTION 3.3 #############"
Rscript scripts\Heaton.R %quick%

echo

echo "All R scripts have been finished!"
pause
exit
