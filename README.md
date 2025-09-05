# CFU_MLE_estimator

A simple scripts in Python and R implementing the approach [used here](https://journals.asm.org/doi/10.1128/spectrum.03946-23))
for estimating CFUs from a dilution series. The authors provide an online calculator [here](https://huggingface.co/spaces/KMichaelMartini/CFUestimator) but 
I wanted to code the approach myself to make sure I was understanding it correctly. I coded it in R and also Python. The author's original python code is 
available here: https://github.com/KMichaelMartini/CFUestimators/tree/main

## Use

To use the scripts you should first count bacterial colonies! Ideally you do this at multiple dilutions but one dilution also works. You should also aim to 
measure the maximum number of technical replicates possible for each dilution. We typically do 3 replicates for a 4 different dilutions in a 96 well format on 
a Nunc Omnitray ([see protocol here](https://github.com/slhogle/OT2_plate_spotting)). You can also perform two replicates at 6 different dilutions using the same 96-well format.

You record these colony counts in a file formatted like `/_data_raw/CFU_test.csv`. This file contains 3 columns: 

1. Column 1 contains the actual observed colony counts
2. Column 2 contains the dilution at which those colonies were counted
3. Column 3 contains some kind of grouping variable (e.g., sample identifier)

There is no need to distinguish replicates just include them in separate rows. When you run the script it will output a CFU/ml estimate based on a truncated Poisson method and a
most probable number approach. The parameters you need to edit include the volume of sample plated or spotted (this volume needs to be the same for the
whole dilution series) and different cutoffs. The crowding cutoff `N` for the censored Poisson approach should be a count value less than you suspect crowding
might start to impact the results (e.g., colonies merging together). For the 2 ul spotting assay a good value for this is 40 to 50. For a standard 100 mm 
Petri dish a good value is 300. The crowding cutoff `NMLE` sets a crowding threshold for the MPN estimator. A good starting value is the ratio of the total 
area of the plate/spot to the average colony size. For a 100 mm petri  dish a good value is 5000. For a 2 ul spot, a decent value to try is 100. The volume 
`V` should be in ml so that the reported CFU density is in units of CFU/ml.

## Availability

Data and code in this GitHub repository (<https://github.com/slhogle/CFU_MLE_estimator>) are provided under [GNU AGPL3](https://www.gnu.org/licenses/agpl-3.0.html).

## Reproducibility

The project uses [`renv`](https://rstudio.github.io/renv/index.html) to create a reproducible environment to execute the code in this project. [See here](https://rstudio.github.io/renv/articles/renv.html#collaboration) for a brief overview on collaboration and reproduction of the entire project. 
To get up and running from an established repository, you could do:

``` r
install.packages("renv")
renv::restore()
```

To initiate `renv` for a new project:

``` r
# if on linux set cran here to download binaries
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))
install.packages("renv")
# initialize
renv::init()
# install some new packages
renv::install("tidyverse")
# record those packages in the lockfile
renv::snapshot()
```
