all: INSTALL_DEPENDS DATA Poisson_sim Negbinom_sim Heaton Sydney Am Chicago MODIS

clean:
	rm -f intermediates/*
	rm -f results/*

# A few targets that are not intended for public use 
# (i.e., they are intended to be run by the original author at submission time)
FIND_DEPENDS: scripts/Dependencies_find.R
	Rscript scripts/Dependencies_find.R

HTML: 
	# Remove all.html if it exists (avoid duplication of entries)
	rm -f all.html
	touch all.html
	# Make the shell script executable:
	chmod +x writehtml.sh
	# populate all.html with images and tables
	./writehtml.sh	

# Targets that are intended to be used by everyone
INSTALL_DEPENDS: scripts/Dependencies_install.R
	Rscript scripts/Dependencies_install.R


DATA: DATA.sh
	chmod +x DATA.sh
	if [ ! -d "data/Sydney_shapefiles" ]; then mkdir data/Sydney_shapefiles; fi
	./DATA.sh


# Targets that do NOT need a "quick" version:
Poisson_sim_OBJECTS = results/Poisson_sim.png results/Poisson_sim_true_process_and_data.png results/Poisson_nres_comparison.csv
Poisson_sim: $(Poisson_sim_OBJECTS) 
$(Poisson_sim_OBJECTS): scripts/Poisson_sim.R
	Rscript scripts/Poisson_sim.R


Negbinom_sim_OBJECTS = results/Negbinom_sim_data.png results/Negbinom_sim_BAU_predictions.png results/Negbinom_sim_arbitrary_polygon_predictions.png results/Negbinom_sim.csv
Negbinom_sim: $(Negbinom_sim_OBJECTS)
$(Negbinom_sim_OBJECTS): scripts/Negbinom_sim.R 
	Rscript scripts/Negbinom_sim.R
scripts/Negbinom_sim.R: scripts/Negbinom_SpatialPolygon_fns.R 


Am: results/Am_comparison.png
results/Am_comparison.png: scripts/Am.R 
	Rscript scripts/Am.R


Sydney_OBJECTS = results/Sydney_training_data.png results/Sydney_SA1_predictions.png results/Sydney_SA3_predictions.png results/Sydney_SA1_coverage.csv
Sydney: $(Sydney_OBJECTS)
$(Sydney_OBJECTS): scripts/Sydney.R 
	Rscript scripts/Sydney.R
scripts/Sydney.R: scripts/Sydney_prep.R


# Targets that need a "quick" version. These targets present some challenges. 
# If we first run "Heaton_quick"", it will produce Heaton_FRKv2.csv. Then, when 
# we try to run "Heaton", it will see that Heaton_FRKv2.csv exists and was 
# created after scripts/Heaton.R: Hence, make will think that nothing is to be
# done and it simply wont re-run the code. A way to circumvent this problem is 
# to use the -B switch to make (long form --always-make), which tells make to 
# disregard timestamps and make the specified targets. This may defeat the 
# purpose of using make, but it's what we need for this submission. 
# Another solution is to use "phony" targets: See 
# https://www.gnu.org/software/make/manual/html_node/Phony-Targets.html. 
# An even simpler workaround is to "touch" a dependency at the end of a run. 
# This is the approach that we take. 

Heaton: results/Heaton_FRKv2.csv 
results/Heaton_FRKv2.csv: scripts/Heaton.R
	Rscript scripts/Heaton.R $(quick)
	if [ "$(quick)" == "TRUE" ]; then touch scripts/Heaton.R; fi


MODIS_OBJECTS = results/MODIS_data.png results/MODIS_MAR_predictions.png results/MODIS_block_predictions.png results/MODIS_ROC.png
MODIS: $(MODIS_OBJECTS)
$(MODIS_OBJECTS): scripts/MODIS.R 
	Rscript scripts/MODIS.R $(quick)
	if [ $(quick) = TRUE ]; then touch scripts/MODIS.R; fi
scripts/MODIS.R: scripts/MODIS_analysis.R 
scripts/MODIS_analysis.R: scripts/MODIS_modelling_fns


Chicago_OBJECTS = results/Chicago_data_pred_uncertainty.png results/Chicago_focused_CAs_time_series.png results/Chicago_focused_CAs_predictive_distributions.png results/Chicago_coverage_and_MAPE.csv
Chicago: $(Chicago_OBJECTS)
$(Chicago_OBJECTS): scripts/Chicago.R
	Rscript scripts/Chicago.R $(quick)
	if [ "$(quick)" == "TRUE" ]; then touch scripts/Chicago.R; fi
scripts/Chicago.R: scripts/Chicago_prep.R