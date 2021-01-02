## Supporting code for "A model for the fragmentation kinetics of crumpled thin sheets"
Jovana Andrejevic, Lisa M. Lee, Shmuel M. Rubinstein, and Chris H. Rycroft

[![DOI](https://zenodo.org/badge/314614182.svg)](https://zenodo.org/badge/latestdoi/314614182)

#### Requirements
Please refer to ```requirements.txt``` for the Python version and package versions used.

Linking to [LAPACK](http://www.netlib.org/lapack/) libraries is needed for running the simulation test written in C++.

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

#### Analysis
- ```plot_dist.py``` - Facet area distributions.
- ```plot_len.py``` - Segment length distributions.
- ```plot_shape.py``` - Shape parameter function.
- ```plot_params.py``` - Various derived parameters from analytical distribution.
- ```plot_scaling.py``` - Scaling results comparison between empirical and derived models.
- ```plot_rates.py``` - Fragmentation rates.
