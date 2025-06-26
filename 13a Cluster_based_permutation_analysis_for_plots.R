library(tidyverse)
library(mc2d) # for pempiricalD()
# For parallelization
library(parallel)
library(doSNOW)


# paths
path <- paste0(getwd(),"/")

#load("rdata/labexpvr_physio_processed.RData")

# Functions --------------------------------------------------------------------

## General functions -----------------------------------------------------------

# Find clusters in the vector of t-statistics. `x` must be a call to
# pt_ttest_statistic_between() or pt_ttest_statistic_within() and `critical_t`
# must be a call to pt_critical_t_between() or pt_critical_t_between(),
# respectively.
# 
# Returns a tibble with the following columns:
# `start`  index of the start of the cluster(s)
# `end`    index of the end of the cluster(s)
# `length` length of the cluster(s)
# (all indices with respect to the series of time points)
clusters_over_time <- function(x, critical_t)
{
  # Find clusters in the t-statistic vector that are above/below the threshold
  x <- rle(abs(x$value) > critical_t)
  
  # Extract only clusters that are above the threshold
  cluster_lengths <- x$lengths[x$values]
  
  # Find the start and end of clusters via the cumsum of cluster lengths. For
  # cluster start, we'll put a 0/FALSE in front, so that we don't run into
  # indexing vec[0]
  cluster_start <- c(0L, cumsum(x$lengths))[which(c(FALSE, x$values)) - 1] + 1L
  cluster_end <- cumsum(x$lengths)[x$values]
  
  tibble(cluster = seq_along(cluster_lengths), start = cluster_start,
         end = cluster_end, length = cluster_lengths)
}

# clusters <- clusters_over_time(t, pt_critical_t_between(tmp, id = "id"))

## Functions for within-subject comparisons (two conditions) -------------------

# This is an implementation of a cluster-based permutation test for within-
# subject comparisons (two conditions) of 2d data (signal x time)
# 
# The data needs to be in a tibble in long format with the following columns and
# must not have any missing data (NAs or missing rows, i.e., each subject needs
# to have 2 * "number of sample" rows):
# 
# (1) A column containing the dependent variable (dv) (e.g., heart rate, skin
#     conductance)
# (2) A column indicating the time points / samples (time)
# (3) A column containing the within-group identifier (within)
# (4) A column containing the participant identifier (id)
#
# All functions expect the column names as character.

# Function for calculating the t-statistic of a within-subject t-test. Used
# instead of `t.test(...)$statistic` for performance increase
t_statistic_within <- function(x, y)
{
  x <- x - y
  
  mean(x) / sqrt(var(x) / length(x))
}

# Calculate the t-statistic (within-subject comparison) for each time point
#data = eda_within; dv = "value"; within = "cue"; time = "bin"
pt_ttest_statistic_within <- function(data, dv, within, time)
{
  data <- data |>
    # Convert to wide format (i.e., values of the `dv` of both within-subject
    # conditions are stored in separate columns
    pivot_wider(names_from = !!within, values_from = !!dv) |>
    # For each sample in the signal ...
    group_by(!!sym(time)) |>
    # ... calculate the t-statistic
    summarize(value = t_statistic_within(
      # From the current subset of the tibble (i.e., per sample), extract the
      # column containing data of the first ...
      .data[[as.character(unique(data[[within]])[1])]],
      # ... and second level of the within-subjects condition
      .data[[as.character(unique(data[[within]])[2])]])) |>
    ungroup()
}

# Calculate the critical t-value for a paired t-test for the given sample size
pt_critical_t_within <- function(data, id, alpha = .05)
{
  qt(1 - alpha / 2, length(unique(data[[id]])) - 1)
}

# Calculation of the null distribution (H0) of cluster lengths for the
# within-subjects comparison
#data = tmp;  dv="scl" ; within = "condition";  time = "timebin"; id = "ID"; nperm = 100
#data = tmp; dv = "value"; within = "cue"; time = "time"; id = "ID"; nperm = 100
#tmp <- alpha_within[[j]]
#tmp <- fnum_within[[j]]
pt_null_distribution_within <- function(data, dv, within, time, id, nperm = 100)
{
  # Create a vector to store the null distribution of cluster lengths
  null_distribution <- vector("integer", nperm)
  
  critical_t <- pt_critical_t_within(data, id = id)
  
  for (i in seq_len(nperm))
  {
    # Make sure that `data` is arranged correctly
    data <- data |> arrange(!!!syms(c(id, within, time)))
    
    # Create a random permutation of the data (= permute the `within` column)
    data <-
      data |>
      # For each participant ...
      group_by(!!sym(id)) |>
      # ... shuffle the within-subject condition and repeat each for the number
      # of samples in the signal
      mutate("{within}" := rep(sample(unique(data[[within]])),
                               each = length(unique(data[[time]])))) |>
      ungroup()
    
    # Find clusters in the random permutation
    t <- pt_ttest_statistic_within(data, dv = dv, within = within, time = time)
    clusters <- clusters_over_time(t, critical_t)
    
    # If clusters were found, take the one with the maximum length and save it
    null_distribution[i] <- if (nrow(clusters) > 0) max(clusters$length) else 0
  }
  
  null_distribution
}

# Prepare data -----------------------------------------------------------------

# For between-subject comparisons, we prepare the data on the fly

# For within-subject comparisons, prepare data frames for comparisons between
# specific phases of the experiment

eda_long <-  readRDS("EDA_long.RData") %>%
  rename(bin = timebin, value = scl, cue = condition)
hr <- read_rds("heart_df.rds")%>%
  filter(problem==0)%>%
  rename(bin = time)%>%
  group_by(vp,bin,cue, ActPass, Intensity)%>%
  summarize(
    value = mean(HRchange, na.rm = T)
  )%>% ungroup()
#%>%
 # filter(as.numeric(bin) < 11)
pupil_long <-  read.csv2(paste0(path.pupil, "pupil_long.csv"))%>% select(-1)%>%
  mutate(cue = as.character(cue))
cb <- read.csv2(paste0(path, "Data/Tobii/CB_Eye_Data_Mean.csv"))
fdur <- read.csv2(paste0(path, "Data/Tobii/FixDur_Eye_Data_Mean.csv"))
fnum <- read.csv2(paste0(path, "Data/Tobii/FixNum_Eye_Data_Mean.csv"))
alpha_supression <- read.csv2("alpha_long.csv") %>%
  mutate(cue = case_when(ActPass == "Active" & Intensity == "High" ~ "HighThreatActive",
                         ActPass == "Active" & Intensity == "Low" ~ "LowThreatActive",
                         ActPass == "Passive" & Intensity == "High" ~ "HighThreatPassive",
                         ActPass == "Passive" & Intensity == "Low" ~ "LowThreatPassive"))%>%
  arrange(ID, cue)%>%
  select(ID,cue,time,value, ActPass, Intensity)%>%
  mutate(bin = time/2)


eda_within <- list(
  
  # Passive: High Intensity vs. Low Intensity
  eda.passive.high_vs_low = eda_long |> filter(ActPass %in% c("Passive")) |>
    select(-ActPass, -Intensity),
  
  # Active: High Intensity vs. Low Intensity
  eda.active.high_vs_low = eda_long |> filter(ActPass %in% c("Active"))|>
    select(-ActPass, -Intensity),
  
  # High Intensity: Active vs. Passive
  eda.high.act_vs_pass = eda_long |> filter(Intensity %in% c("High"))|>
    select(-ActPass, -Intensity),
  
  # Low Intensity: Active vs. Passive
  eda.low.act_vs_pass = eda_long |> filter(Intensity %in% c("Low"))|>
    select(-ActPass, -Intensity)
)


hr_within <- list(

  # Passive: High Intensity vs. Low Intensity
  hr.passive.high_vs_low = hr |> filter(ActPass %in% c("Passive")) |>
    select(-ActPass, -Intensity),
  
  # Active: High Intensity vs. Low Intensity
  hr.active.high_vs_low = hr |> filter(ActPass %in% c("Active"))|>
    select(-ActPass, -Intensity),
  
  # High Intensity: Active vs. Passive
  hr.high.act_vs_pass = hr |> filter(Intensity %in% c("High"))|>
    select(-ActPass, -Intensity),
  
  # Low Intensity: Active vs. Passive
  hr.low.act_vs_pass = hr |> filter(Intensity %in% c("Low"))|>
    select(-ActPass, -Intensity)
)


pupil_within <- list(

  # Passive: High Intensity vs. Low Intensity
  pupil.passive.high_vs_low = pupil_long |> filter(ActPass %in% c("Passive")) |>
    select(-ActPass, -Intensity),

  # Active: High Intensity vs. Low Intensity
  pupil.active.high_vs_low = pupil_long |> filter(ActPass %in% c("Active"))|>
    select(-ActPass, -Intensity),

  # High Intensity: Active vs. Passive
  pupil.high.act_vs_pass = pupil_long |> filter(Intensity %in% c("High"))|>
    select(-ActPass, -Intensity),

  # Low Intensity: Active vs. Passive
  pupil.low.act_vs_pass = pupil_long |> filter(Intensity %in% c("Low"))|>
    select(-ActPass, -Intensity)

)

cb_within <- list(
  
  # Passive: High Intensity vs. Low Intensity
  cb.passive.high_vs_low = cb |> filter(ActPass %in% c("Passive")) |>
    select(-ActPass, -Threat),
  
  # Active: High Intensity vs. Low Intensity
  cb.active.high_vs_low = cb |> filter(ActPass %in% c("Active"))|>
    select(-ActPass, -Threat),
  
  # High Intensity: Active vs. Passive
  cb.high.act_vs_pass = cb |> filter(Threat %in% c("High"))|>
    select(-ActPass, -Threat),
  
  # Low Intensity: Active vs. Passive
  cb.low.act_vs_pass = cb |> filter(Threat %in% c("Low"))|>
    select(-ActPass, -Threat)
  
)

fdur_within <- list(
  
  # Passive: High Intensity vs. Low Intensity
  fdur.passive.high_vs_low = fdur |> filter(ActPass %in% c("Passive")) |>
    select(-ActPass, -Threat),
  
  # Active: High Intensity vs. Low Intensity
  fdur.active.high_vs_low = fdur |> filter(ActPass %in% c("Active"))|>
    select(-ActPass, -Threat),
  
  # High Intensity: Active vs. Passive
  fdur.high.act_vs_pass = fdur |> filter(Threat %in% c("High"))|>
    select(-ActPass, -Threat),
  
  # Low Intensity: Active vs. Passive
  fdur.low.act_vs_pass = fdur |> filter(Threat %in% c("Low"))|>
    select(-ActPass, -Threat)
  
)

fnum_within <- list(
  
  # Passive: High Intensity vs. Low Intensity
  fnum.passive.high_vs_low = fnum |> filter(ActPass %in% c("Passive")) |>
    select(-ActPass, -Threat),
  
  # Active: High Intensity vs. Low Intensity
  fnum.active.high_vs_low = fnum |> filter(ActPass %in% c("Active"))|>
    select(-ActPass, -Threat),
  
  # High Intensity: Active vs. Passive
  fnum.high.act_vs_pass = fnum |> filter(Threat %in% c("High"))|>
    select(-ActPass, -Threat),
  
  # Low Intensity: Active vs. Passive
  fnum.low.act_vs_pass = fnum |> filter(Threat %in% c("Low"))|>
    select(-ActPass, -Threat)
  
)

alpha_within <- list(
  
  # Passive: High Intensity vs. Low Intensity
  alpha.passive.high_vs_low = alpha_supression |> filter(ActPass %in% c("Passive")) |>
    select(-ActPass, -Intensity, -time),
  
  # Active: High Intensity vs. Low Intensity
  alpha.active.high_vs_low = alpha_supression |> filter(ActPass %in% c("Active"))|>
    select(-ActPass, -Intensity, -time),
  
  # High Intensity: Active vs. Passive
  alpha.high.act_vs_pass = alpha_supression |> filter(Intensity %in% c("High"))|>
    select(-ActPass, -Intensity, -time),
  
  # Low Intensity: Active vs. Passive
  alpha.low.act_vs_pass = alpha_supression |> filter(Intensity %in% c("Low"))|>
    select(-ActPass, -Intensity, -time)
  
)

# Calculate null distributions -------------------------------------------------
iterations = 1000

# Use parallelization to run the permutations. As this was written for Windows,
# the doSNOW backend is used in the current implementation.

## Prepare parallelization ----------------------------------------------------

cl <- makeCluster(detectCores() - 1)
registerDoSNOW(cl)

## EDA within ------------------------------------------------------------------

eda_null_dist_within <- setNames(vector("list", length(eda_within)),
                                 names(eda_within))

for (j in seq_along(eda_within)){
  # Extract the data of the current comparison
  tmp <- eda_within[[j]] 
  
  # Calculate null distribution
  x1 <- foreach(i = 1:length(cl), .combine = "c", .packages = "tidyverse") %dopar% {
    y <- pt_null_distribution_within(tmp, dv = "value", within = "cue",
                                     time = "bin", id = "vp", nperm = ceiling(iterations/length(cl)))
    return(y)
  }
  
  eda_null_dist_within[[j]] <- x1[1:1000]
}

# Create a tibble with the 95% quantile of cluster lengths under H0 per phase
eda_critical_length_null_within <-
  tibble(cue = names(eda_null_dist_within),
         critical_length = map_dbl(eda_null_dist_within, quantile, .95))

## HR within -------------------------------------------------------------------

hr_null_dist_within <- setNames(vector("list", length(hr_within)),
                                names(hr_within))

for (j in seq_along(hr_within)){
  # Extract the data of the current comparison
  tmp <- hr_within[[j]]
  
  # Calculate null distribution
  x1 <- foreach(i = 1:length(cl), .combine = "c", .packages = "tidyverse") %dopar% {
    y <- pt_null_distribution_within(tmp, dv = "value", within = "cue",
                                     time = "bin", id = "vp", nperm = ceiling(iterations/length(cl)))
    return(y)
  }
  
  hr_null_dist_within[[j]] <- x1[1:iterations]
}

# Create a tibble with the 95% quantile of cluster lengths under H0 per phase
hr_critical_length_null_within <-
  tibble(cue = names(hr_null_dist_within),
         critical_length = map_dbl(hr_null_dist_within, quantile, .95))


## Pupil within -------------------------------------------------------------------

pupil_null_dist_within <- setNames(vector("list", length(pupil_within)),
                                   names(pupil_within))

for (j in seq_along(pupil_within)){
  # Extract the data of the current comparison
  tmp <- pupil_within[[j]]
  
  # Calculate null distribution
  x1 <- foreach(i = 1:length(cl), .combine = "c", .packages = "tidyverse") %dopar% {
    y <- pt_null_distribution_within(tmp, dv = "value", within = "cue",
                                     time = "bin", id = "code", nperm = ceiling(iterations/length(cl)))
    return(y)
  }
  
  pupil_null_dist_within[[j]] <- x1[1:iterations]
}

# Create a tibble with the 95% quantile of cluster lengths under H0 per phase
pupil_critical_length_null_within <-
  tibble(cue = names(pupil_null_dist_within),
         critical_length = map_dbl(pupil_null_dist_within, quantile, .95))


## CB within -------------------------------------------------------------------

cb_null_dist_within <- setNames(vector("list", length(cb_within)),
                                   names(cb_within))

for (j in seq_along(cb_within)){
  # Extract the data of the current comparison
  tmp <- cb_within[[j]]
  
  # Calculate null distribution
  x1 <- foreach(i = 1:length(cl), .combine = "c", .packages = "tidyverse") %dopar% {
    y <- pt_null_distribution_within(tmp, dv = "value", within = "cue",
                                     time = "bin", id = "vp", nperm = ceiling(iterations/length(cl)))
    return(y)
  }
  
  cb_null_dist_within[[j]] <- x1[1:iterations]
}

# Create a tibble with the 95% quantile of cluster lengths under H0 per phase
cb_critical_length_null_within <-
  tibble(cue = names(cb_null_dist_within),
         critical_length = map_dbl(cb_null_dist_within, quantile, .95))


## FNum within -------------------------------------------------------------------

fnum_null_dist_within <- setNames(vector("list", length(fnum_within)),
                                names(fnum_within))

for (j in seq_along(fnum_within)){
  # Extract the data of the current comparison
  tmp <- fnum_within[[j]]
  
  # Calculate null distribution
  x1 <- foreach(i = 1:length(cl), .combine = "c", .packages = "tidyverse") %dopar% {
    y <- pt_null_distribution_within(tmp, dv = "value", within = "cue",
                                     time = "bin", id = "vp", nperm = ceiling(iterations/length(cl)))
    return(y)
  }
  
  fnum_null_dist_within[[j]] <- x1[1:iterations]
}

# Create a tibble with the 95% quantile of cluster lengths under H0 per phase
fnum_critical_length_null_within <-
  tibble(cue = names(fnum_null_dist_within),
         critical_length = map_dbl(fnum_null_dist_within, quantile, .95))

## FDur within -------------------------------------------------------------------

fdur_null_dist_within <- setNames(vector("list", length(fdur_within)),
                                names(fdur_within))

for (j in seq_along(fdur_within)){
  # Extract the data of the current comparison
  tmp <- fdur_within[[j]]
  
  # Calculate null distribution
  x1 <- foreach(i = 1:length(cl), .combine = "c", .packages = "tidyverse") %dopar% {
    y <- pt_null_distribution_within(tmp, dv = "value", within = "cue",
                                     time = "bin", id = "vp", nperm = ceiling(iterations/length(cl)))
    return(y)
  }
  
  fdur_null_dist_within[[j]] <- x1[1:iterations]
}

# Create a tibble with the 95% quantile of cluster lengths under H0 per phase
fdur_critical_length_null_within <-
  tibble(cue = names(fdur_null_dist_within),
         critical_length = map_dbl(fdur_null_dist_within, quantile, .95))

## Alpha within -------------------------------------------------------------------

alpha_null_dist_within <- setNames(vector("list", length(alpha_within)),
                                  names(alpha_within))

for (j in seq_along(alpha_within)){
  # Extract the data of the current comparison
  tmp <- alpha_within[[j]]
  
  # Calculate null distribution
  x1 <- foreach(i = 1:length(cl), .combine = "c", .packages = "tidyverse") %dopar% {
    y <- pt_null_distribution_within(tmp, dv = "value", within = "cue",
                                     time = "bin", id = "ID", nperm = ceiling(iterations/length(cl)))
    return(y)
  }
  
  alpha_null_dist_within[[j]] <- x1[1:iterations]
}

# Create a tibble with the 95% quantile of cluster lengths under H0 per phase
alpha_critical_length_null_within <-
  tibble(cue = names(alpha_null_dist_within),
         critical_length = map_dbl(alpha_null_dist_within, quantile, .95))

### Save distributions and stop cluster ----------------------------------------

save(hr_null_dist_within, hr_critical_length_null_within,
     eda_null_dist_within, eda_critical_length_null_within,
     pupil_null_dist_within, pupil_critical_length_null_within,
     cb_null_dist_within, cb_critical_length_null_within,
     fnum_null_dist_within, fnum_critical_length_null_within,
     fdur_null_dist_within, fdur_critical_length_null_within,
     alpha_null_dist_within, alpha_critical_length_null_within,
     file = "physio_pt_null_distributions.RData")

stopCluster(cl)

## Cluster-based permutation tests ---------------------------------------------

load("physio_pt_null_distributions.RData")

# Within-group comparison (compare different phases) ---------------------------

## Cluster-based permutation tests ---------------------------------------------

## Cluster-based permutation tests (within) ------------------------------------

### EDA ------------------------------------------------------------------------

# Find clusters of differences in the EDA signal in each phase of the experiment
eda_clusters_within <-
  eda_within |>
  map(pt_ttest_statistic_within, dv = "value", within = "cue", time = "bin") |>
  map(clusters_over_time, critical_t = pt_critical_t_within(eda_long, id = "vp")) |>
  # Bind to a single tibble
  bind_rows(.id = "cue") |>
  # Add critical cluster length under H0 (95% quantile)
  left_join(eda_critical_length_null_within, by = "cue") |>
  mutate(start = as.numeric(start),
         end = as.numeric(end)) %>%
  # Calculate the p-value for each cluster
  rowwise() |>
  mutate(p = 1 - pempiricalD(length, eda_null_dist_within[[cue]])) |>
  ungroup()

### HR -------------------------------------------------------------------------

# Find clusters of differences in the EDA signal in each phase of the experiment
hr_clusters_within <-
  hr_within |>
  map(pt_ttest_statistic_within, dv = "value", within = "cue", time = "bin") |>
  map(clusters_over_time, critical_t = pt_critical_t_within(hr, id = "vp")) |>
  # Bind to a single tibble
  bind_rows(.id = "cue") |>
  # Add critical cluster length under H0 (95% quantile)
  left_join(hr_critical_length_null_within, by = "cue") |>
  # Calculate the p-value for each cluster
  rowwise() |>
  mutate(p = 1 - pempiricalD(length, hr_null_dist_within[[cue]])) |>
  ungroup()

### Pupil -------------------------------------------------------------------------

# Find clusters of differences in the EDA signal in each phase of the experiment
pupil_clusters_within <-
  pupil_within |>
  map(pt_ttest_statistic_within, dv = "value", within = "cue", time = "bin") |>
  map(clusters_over_time, critical_t = pt_critical_t_within(pupil_long, id = "code")) |>
  # Bind to a single tibble
  bind_rows(.id = "cue") |>
  # Add critical cluster length under H0 (95% quantile)
  left_join(pupil_critical_length_null_within, by = "cue") |>
  # Calculate the p-value for each cluster
  rowwise() |>
  mutate(p = 1 - pempiricalD(length, pupil_null_dist_within[[cue]])) |>
  mutate(start = start/2)|>
  mutate(end = end/2)|>
  ungroup()

### cb -------------------------------------------------------------------------

# Find clusters of differences in the EDA signal in each phase of the experiment
cb_clusters_within <-
  cb_within |>
  map(pt_ttest_statistic_within, dv = "value", within = "cue", time = "bin") |>
  map(clusters_over_time, critical_t = pt_critical_t_within(cb, id = "vp")) |>
  # Bind to a single tibble
  bind_rows(.id = "cue") |>
  # Add critical cluster length under H0 (95% quantile)
  left_join(cb_critical_length_null_within, by = "cue") |>
  # Calculate the p-value for each cluster
  rowwise() |>
  mutate(p = 1 - pempiricalD(length, cb_null_dist_within[[cue]])) |>
  mutate(start = start+2)|>
  mutate(end = end+2)|>
  ungroup()

### fnum -------------------------------------------------------------------------

# Find clusters of differences in the EDA signal in each phase of the experiment
fnum_clusters_within <-
  fnum_within |>
  map(pt_ttest_statistic_within, dv = "value", within = "cue", time = "bin") |>
  map(clusters_over_time, critical_t = pt_critical_t_within(fnum, id = "vp")) |>
  # Bind to a single tibble
  bind_rows(.id = "cue") |>
  # Add critical cluster length under H0 (95% quantile)
  left_join(fnum_critical_length_null_within, by = "cue") |>
  # Calculate the p-value for each cluster
  rowwise() |>
  mutate(p = 1 - pempiricalD(length, fnum_null_dist_within[[cue]])) |>
  mutate(start = start+2)|>
  mutate(end = end+2)|>
  ungroup()


### fdur -------------------------------------------------------------------------

# Find clusters of differences in the EDA signal in each phase of the experiment
fdur_clusters_within <-
  fdur_within |>
  map(pt_ttest_statistic_within, dv = "value", within = "cue", time = "bin") |>
  map(clusters_over_time, critical_t = pt_critical_t_within(fdur, id = "vp")) |>
  # Bind to a single tibble
  bind_rows(.id = "cue") |>
  # Add critical cluster length under H0 (95% quantile)
  left_join(fdur_critical_length_null_within, by = "cue") |>
  # Calculate the p-value for each cluster
  rowwise() |>
  mutate(p = 1 - pempiricalD(length, fdur_null_dist_within[[cue]])) |>
  mutate(start = start+2)|>
  mutate(end = end+2)|>
  ungroup()

### alpha ------------------------------------------------------------------------

# Find clusters of differences in the EDA signal in each phase of the experiment
alpha_clusters_within <-
  alpha_within |>
  map(pt_ttest_statistic_within, dv = "value", within = "cue", time = "bin") |>
  map(clusters_over_time, critical_t = pt_critical_t_within(alpha_supression, id = "ID")) |>
  # Bind to a single tibble
  bind_rows(.id = "cue") |>
  # Add critical cluster length under H0 (95% quantile)
  left_join(alpha_critical_length_null_within, by = "cue") |>
  # Calculate the p-value for each cluster
  rowwise() |>
  mutate(p = 1 - pempiricalD(length, alpha_null_dist_within[[cue]])) |>
  mutate(start = start/2)|>
  mutate(end = end/2)|>
  ungroup()

