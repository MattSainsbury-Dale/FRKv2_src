all: FIND_DEPENDS INSTALL_DEPENDS DATA Poisson_sim Negbinom_sim Heaton MODIS Sydney Am Chicago

FIND_DEPENDS: src/Dependencies_find.R
	Rscript src/Dependencies_find.R

INSTALL_DEPENDS: src/Dependencies_install.R
	Rscript src/Dependencies_install.R

DATA: DATA.sh
	chmod +x DATA.sh
	./DATA.sh

Poisson_sim_OBJECTS = img/Poisson_sim.png img/Poisson_sim_true_process_and_data.png results/Poisson_nres_comparison.csv

Poisson_sim: $(Poisson_sim_OBJECTS) 

$(Poisson_sim_OBJECTS): src/Poisson_sim.R
	Rscript src/Poisson_sim.R

Negbinom_sim_OBJECTS = img/Negbinom_sim_data.png img/Negbinom_sim_BAU_predictions.png img/Negbinom_sim_arbitrary_polygon_predictions.png results/Negbinom_sim.csv

Negbinom_sim: $(Negbinom_sim_OBJECTS)

$(Negbinom_sim_OBJECTS): src/Negbinom_sim.R 
	Rscript src/Negbinom_sim.R

src/Negbinom_sim.R: src/Negbinom_SpatialPolygon_fns.R 

Heaton: results/Heaton_FRKv2.csv

results/Heaton_FRKv2.csv: src/Heaton.R 
	Rscript src/Heaton.R


MODIS_OBJECTS = img/MODIS_data.png img/MODIS_MAR_predictions.png img/MODIS_block_predictions.png img/MODIS_ROC.png

MODIS: $(MODIS_OBJECTS)

$(MODIS_OBJECTS): src/MODIS_control.R 
	Rscript src/MODIS_control.R
	
src/MODIS_control.R: src/MODIS.R 

src/MODIS.R: src/MODIS_modelling_fns


Am: img/Am_data_and_blocks.png img/Am_comparison.png

img/Am_comparison.png: src/Am_comparison_plot.R intermediates/Am_FRK.csv intermediates/Am_georob.csv
	Rscript src/Am_comparison_plot.R

intermediates/Am_FRK.csv: src/Am_FRK.R intermediates/Am_data.csv intermediates/Am_BAUs.rds intermediates/Am_blocks.rds
	Rscript src/Am_FRK.R

intermediates/Am_georob.csv: src/Am_georob.R intermediates/Am_data.csv intermediates/Am_BAUs.rds intermediates/Am_blocks.rds
	Rscript src/Am_georob.R
	
img/Am_data_and_blocks.png: src/Am_data_blocks_plot.R intermediates/Am_blocks.rds intermediates/Am_GZ.csv intermediates/Am_data.csv
	Rscript src/Am_data_blocks_plot.R

intermediates/Am_BAUs.rds: src/Am_BAUs.R intermediates/Am_data.csv intermediates/Am_GZ.csv
	Rscript src/Am_BAUs.R
	
intermediates/Am_blocks.rds: src/Am_blocks.R intermediates/Am_GZ.csv
	Rscript src/Am_blocks.R

intermediates/Am_data.csv: src/Am_data.R data/Am_data.csv intermediates/Am_GZ.csv
	Rscript src/Am_data.R

intermediates/Am_GZ.csv: src/Am_GZ.R
	Rscript src/Am_GZ.R


Sydney_OBJECTS = img/Sydney_training_data.png img/Sydney_SA1_predictions.png img/Sydney_SA3_predictions.png results/Sydney_SA1_coverage.csv

Sydney: $(Sydney_OBJECTS)

$(Sydney_OBJECTS): src/Sydney.R 
	Rscript src/Sydney.R
	
src/Sydney.R: src/Sydney_prep.R



Chicago_OBJECTS = img/Chicago_data_pred_uncertainty.png img/Chicago_focused_CAs_time_series.png img/Chicago_focused_CAs_predictive_distributions.png results/Chicago_coverage_and_MAPE.csv

Chicago: $(Chicago_OBJECTS)

$(Chicago_OBJECTS): src/Chicago.R
	Rscript src/Chicago.R
	
src/Chicago.R: src/Chicago_prep.R

clean:
	rm -f img/*
	rm -f intermediates/*
	rm -f results/*