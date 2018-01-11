# allanhills
Scripts and utilities for assessing the potential for a million-year-old ice-core record from the Allan Hills Blue Ice Area, East Antarctica, using ice-penetrating radar observations, age constraints, and a simple flowline model.

Code developed as part of my dissertation at the University of Washington.

## Dependencies
- Matlab toolboxes, including [Antarctic Mapping Tools](https://www.mathworks.com/matlabcentral/fileexchange/47638-antarctic-mapping-tools), [Cubehelix](https://www.mathworks.com/matlabcentral/fileexchange/43700-cubehelix-colormaps--beautiful--distinct--versatile-), [BedMap2](https://www.mathworks.com/matlabcentral/fileexchange/42353-bedmap2-toolbox-for-matlab), [Antarctic Drainage Basins](https://www.mathworks.com/matlabcentral/fileexchange/47639-antarctic-drainage-basins), [Antarctic boundaries, grounding line, and masks from InSAR](https://www.mathworks.com/matlabcentral/fileexchange/60246-antarctic-boundaries--grounding-line--and-masks-from-insar), and [export_fig](https://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig)

## To run code
- `setup_paths.m` - sets up the necessary paths to run scripts (some paths will need to be adapted for the user's environment)
- `scripts/make_plots.m` - loads all data and must be run before other scripts

## allanhills/scripts
- `grinsted_mode.m` - runs a model developed by [Grinsted et al [2003]](http://onlinelibrary.wiley.com/doi/10.1029/2003GL017957/full) for various combinations of snow accumulation rate and surface velocity decreases during glacial periods, calculates the misfit between the model and the dated isochrones, and plots the results
- `make_plots.m` - loads, processes, and plots the radar observations, age constraints, and surface velocities 
- `internal_stratigraphy.m` - calculates the height above the bed for different dated isochrones and creates a 3D plot to look at the 3D internal stratigraphy
- `layer_differences.m` - calculates layer thickness anomalies along a radar profile
- `load_picked_layers.m` - example code for loading picked layers from spicker using `load_splayer.m`
- `plot_potential_core.m` - plots best-fit models for a potential core site and determines the height of 500-ka or 1-Ma ice above the bed

## allanhills/util
- some of these scripts were authored by others

## allanhills/ages and allanhills/velocities
- data files for age constraints from vertical and horizontal ice cores from [Spaulding et al [2013]](http://www.sciencedirect.com/science/article/pii/S003358941300080X) and for stake measurements from [Spaulding et al [2012]](https://www.igsoc.org/journal/58/208/j11j176.pdf)

## allanhills/pickedlayers
- picked radar layers, which are loaded in other scripts
