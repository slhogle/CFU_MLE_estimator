library(tidyverse)
library(here)

# Global vars -------------------------------------------------------------

# set the poisson threshold N. This should be a value over which you don't trust the counts
# for the 2 ul spotting assay a good value for this is 50
N <- 50

# set the MLE threshold NMLE. A good starting estimate is the ratio of the total area of 
# the plate/spot to the average colony size. Used for calculating MPN. For a 100 mm petri 
# dish a good value is 5000. For a 2 ul spot, a decent value is 100
NMLE <- 100

# the output will be expressed in terms of absolute CFU in the amount that was plated or 
# spotted. It is very important that the same spotting/plating volume be used for
# all dilutions. To convert to CFU/ml divide R by the spot/plating volume. E.g., if these
# are counts from the 96-well spot assay the value is usually 0.002 or 0.0025
sampling_vol <- 0.002

# Functions ---------------------------------------------------------------

find_poisson_cutoff <- function(counts, dil, N, V=1){
  mask <- counts<N
  r_p <- sum(counts[mask])/(sum(dil[mask])*V)
  r_p_std <- sqrt(r_p*r_p/sum(counts[mask]))
  tibble(r = r_p, stderr = r_p_std, type = "poisson_cutoff")
}

find_MLE <- function(counts, dil, N, V=1){
  f <- function(x) {
    sum(counts*dil/(N*(1-exp(-x*dil*V/N))))-sum(dil)
  }
  if (any(N<counts)) {
    tibble(r = Inf, stderr = Inf, type = "ML_estimator")
  } else {
    root <- uniroot(f, interval = c(1, 1e9))
    r_mle <- root$root
    p0 <- exp(-r_mle*dil*V/N)
    invVar <- sum(dil*dil*V*V*counts*p0/(N*N*(1-p0)**2))
    estVar <- 1/invVar
    r_mle_std <- sqrt(estVar)
    tibble(r = r_mle, stderr = r_mle_std, type = "ML_estimator")
  }
}

find_estimators <- function(counts, dil, N, NMLE) {
  bind_rows(find_poisson_cutoff(counts, dil, N),
            find_MLE(counts, dil, NMLE))
}

# Run ---------------------------------------------------------------------

cfu_df <- readr::read_csv(here::here("_data_raw", "CFU_test.csv"), col_names = c("count", "dilution", "group"))

cfu_df %>% 
  nest(data = c(count, dilution)) %>% 
  mutate(r = map(data, \(x) find_estimators(x$count, x$dilution, N=N, NMLE=NMLE))) %>% 
  unnest(r) %>%
  mutate(CFU_ml = r/sampling_vol, 
         CFU_ml_stderr = stderr/sampling_vol) %>% 
  dplyr::select(group, CFU_ml, CFU_ml_stderr, type)

