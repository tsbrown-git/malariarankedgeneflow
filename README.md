# Example code to accompany "Distinguishing gene flow between malaria parasite populations"

# Data

*example_ms.tab*: example command line for generating coalescent simulated sequence data using msprime (https://msprime.readthedocs.io/en/stable/). Generates simulated sequence data for a metapopulations with 10 demes and migration rates, recombination rate, and mutation rate as described in Methods for M=0.1 and Ne=100.

# Code
*R_PRC_simulated.R* : Calculates proportion ranked correctly (PRC) values using R_widehat (a population-based summary of relatedness, based on identity by descent) for all pairs location pairs for one selected replicate from ms-simulated sequence data with 10 demes (which can be generated using *example_ms.tab*). Outputs PRC values for all pairs of location pairs (=990) using a specified number of markers (*p*) for each number of individuals (*n*) in 5,10,...,95,100. Required options: --stem (stem for stem for multipopulation ms file, without .txt extension); --path (path to ms file); --ncores (number of cores to use); --rep (which simulation replicate to use); --snps (number of snps to use). Requires hmmIBD (https://github.com/glipsnort/hmmIBD) and assumes this is installed in working directory. Requires R packages optparse, doParallel, stringr, readr, reshape2, msr (https://github.com/vsbuffalo/msr), readr, and dplyr.

*FST_PRC_simlated.R* : Calculates PRC values using Weir and Cockerhams FST estimator for all pairs location pairs for one selected replicate from ms-simulated sequence data with 10 demes (which can be generated using *example_ms.tab*). Outputs PRC values for all pairs of location pairs (=990) using a specified number of markers (*p*) for each number of individuals (*n*) in 5,10,...,95,100. Required options: --stem (stem for stem for multipopulation ms file, without .txt extension); --path (path to ms file); --ncores (number of cores to use); --rep (which simulation replicate to use); --snps (number of snps to use). Requires R packages optparse, doParallel, stringr, readr, reshape2, msr (https://github.com/vsbuffalo/msr), readr, and dplyr.

