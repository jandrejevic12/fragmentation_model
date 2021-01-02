## Supporting code for "A model for the fragmentation kinetics of crumpled thin sheets"
Jovana Andrejevic, Lisa M. Lee, Shmuel M. Rubinstein, and Chris H. Rycroft

[![DOI](https://zenodo.org/badge/314614182.svg)](https://zenodo.org/badge/latestdoi/314614182)

#### Requirements
Please refer to ```requirements.txt``` for the Python version and package versions used. \
Linking to [LAPACK](http://www.netlib.org/lapack/) library is needed for running the simulation test written in C++.

### Data Availability
Data collected for the study of [*Gottesman et al.*][state_var] was initially processed to delineate individual facets and their properties, and the resulting data structures saved in a separate file for analysis. The post-processed data is provided as ```facet_data.zip``` and may be downloaded to run all analysis scripts in this repository. The raw data from the work of [*Gottesman et al.*][state_var] can be provided from the authors upon reasonable request.

#### Imports
- ```math_imports.py``` - Import special math functions and define custom distribution classes.
- ```plot_imports.py``` - Import plotting libraries and define custom plotting setup.
- ```image_imports.py``` - Import image processing functions.
- ```setup_imports.py``` - Set up directories and common variables.

#### Numerical simulation
- ```frag_test.cc``` - Numerical integration of the fragmentation rate equation.
- ```frag_plot.py``` - Plot simulation results.

#### Processing
- ```facet_utils.py``` - Routines for processing facet segmentation data.
- ```predict_rates.py``` - Predict fragmentation rates from segmentation data.

#### Analysis
- ```plot_dist.py``` - Facet area distributions.
- ```plot_len.py``` - Segment length distributions.
- ```plot_shape.py``` - Shape parameter function.
- ```plot_params.py``` - Various derived parameters from analytical distribution.
- ```plot_scaling.py``` - Scaling results comparison between empirical and derived models.

[state_var]: https://doi.org/10.1038/s42005-018-0072-x "Gottesman, O., Andrejevic, J., Rycroft, C. H., & Rubinstein, S. M. (2018). A state variable for crumpled thin sheets. Communications Physics, 1(1), 1-7. Available from: https://doi.org/10.1038/s42005-018-0072-x."
