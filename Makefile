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