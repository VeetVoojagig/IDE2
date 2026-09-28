## IDE SCRIPT ------------------------------------------------------------------
# By Mark S Gilthorpe
#
# Challenges in evaluating change in heteroscedastic data: a 2x2 simulation
# framework crossing allocation mechanism (random / targeted on baseline)
# against variance structure (homoscedastic / heteroscedastic). 

# Cell A is a single N = 500 draw from the homoscedastic population, every
# child group = 1, with no comparator arm at all (positive control,
# intervention-only test against zero). 
#
# Cell B is TWO independent N = 500 draws from the same heteroscedastic 
# population parametrisation, one per arm, not one population split by coin flip. 
#
# Cells C & D are a draw of two independent arms, each drawn afresh in every
# replication from its own baseline-truncated distribution, thresholded
# IDENTICALLY for treatment & control: the comparator must be drawn from the
# same restricted / truncated population as the treatment arm, since an
# exchangeable comparator has to share the treatment arm's baseline variance
# structure (see "TARGETED-ALLOCATION POOLS" below).
#
# The intervention differential effect (mean shift, & for true_ide, weakened
# tracking - a rescaled raw SD at each follow-up occasion under a fixed
# correlation matrix) is baked directly into the treatment arm's GENERATING
# distribution at draw time (see build_treatment_distribution()) - never
# reconstructed afterwards.
#
# Naming: true_ide in this script corresponds to true_ide in the manuscript (a
# genuine negative IDE); no_ide is unchanged. Likewise Cells A-D here are the
# manuscript's Scenarios A-D, & homoscedastic / heteroscedastic here are the
# manuscript's constant / non-constant outcome variance.
#
# Cell B is the paper's central empirical result (heteroscedastic + random); 
# Cell C is a within-framework confirmatory replication of Beggs et al. - i.e., 
# homoscedastic + targeted - but with much more extreme truncation, which makes
# a vital new point regarding the impacts of extreme truncations; and Cell D is 
# the extension to both complications (heteroscedastic + extremely non-random).
#
# Script structure
#   PART 1  SETUP & GLOBAL PARAMETERS
#   PART 2  FUNCTIONS LIBRARY (definitions only)
#   PART 3  FUNCTIONS VERIFICATION
#   PART 4  SIMULATION, CACHE, & DIAGNOSTICS
#   PART 5  CONFIRMATORY DIAGNOSTIC (Fisher z reference variance)
#   PART 6  FIGURES & TABLES (Type 1 error / Power, no_ide & true_ide)

## #############################################################################
## PART ONE   -- SETUP & GLOBAL PARAMETERS
## PACKAGES   -- LOAD REQUIRED LIBRARIES ---------------------------------------

# Clear the environment; project is named 'IDE'
rm(list = ls()); gc(full = TRUE); closeAllConnections(); Proj_Name <- "IDE"
pkgs     <- c("MASS", "ggplot2", "parallel", "compiler", "psych", "cocor",
              "openxlsx", "grid", "lavaan")
new_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)
for (itn in seq_along(pkgs))
  suppressPackageStartupMessages(library(pkgs[itn], character.only = TRUE))
rm(list = c("itn", "new_pkgs", "pkgs"))

## SETUP      -- RUN-CONTROL FLAGS & SIMULATION SIZES --------------------------

# force_rerun_* : TRUE invalidates the named cache & forces regeneration; lavaan 
# is now the sole estimation route behind MLM. 
# parallel_K    : K cutover above which lapply() switches to parLapply()
# K             : single replication count used by every repeated-replication
#                 pass in the script - the Oldham sanity check (PART THREE),
#                 the MLM pass (PART FOUR), the threshold sweep incl. its MLM
#                 leg (PART FOUR), & the confirmatory diagnostic (PART FIVE).
#                 One value, one place to change it. Set K low (e.g., 100) to
#                 test the full pipeline end-to-end quickly, then raise it for
#                 a real run - this matters most for the MLM sweep leg &
#                 severity plot, since each lavaan fit is expensive & the
#                 "MLM stays near 5%" claim only firms up at a much larger K.
# sweep_K_max   : ceiling on the threshold sweep's replication count (both
#                 legs, PART FOUR): sweep_K = min(K, sweep_K_max). The MLM
#                 sweep leg at K = 1e6 would run for about a week, so it is
#                 held at 1e5; set sweep_K_max <- K to run the sweep at full K.
#                 The PART FIVE diagnostic is Oldham only & still uses full K
# chunk_size    : progress bar, parallel dispatch, & checkpoint granularity
#                 for Oldham/sweep/diagnostic (PART THREE, PART FIVE)
# mlm_chnk_size : same role as chunk_size but for the MLM pass (PART FOUR)
#                 specifically, kept smaller since lavaan growth-model fits
#                 are more expensive per unit than Oldham's closed-form
#                 correlation - a crash mid-run loses at most this many
#                 replications, not up to chunk_size
# sweep_pct     : tail-probability cuts for the sweep (%), evenly spaced
#                 on the log scale; sweep_z is derived automatically

if (!exists("force_rerun_OLDHAM"))  force_rerun_OLDHAM  <- FALSE
if (!exists("force_rerun_MLM"))     force_rerun_MLM     <- FALSE
if (!exists("force_rerun_SWEEP"))   force_rerun_SWEEP   <- FALSE

# skip_appendix : TRUE skips the Appendix pass entirely (alternative
# tracking-preserving true_ide pools, verify_res_alt/mlm_res_alt, Table A1,
# & its xlsx sheet) - this is the single most expensive optional piece of
# the pipeline (the appendix MLM power-only pass alone can run several
# hours at large K) & feeds ONLY Table A1, nothing else - the severity /
# power plot & Table A2 are built entirely from the sweep / diagnostic sections 
# & are completely unaffected by this flag. Set TRUE for a fast run when you 
# only need the plots / Table 2 / Table A2.

if (!exists("skip_appendix"))       skip_appendix       <- FALSE
if (!exists("parallel_K"))          parallel_K          <- 1000
if (!exists("K"))                   K                   <- 1e6
if (!exists("sweep_K_max"))         sweep_K_max         <- 1e5
if (!exists("chunk_size"))          chunk_size          <- 5e3
if (!exists("mlm_chnk_size"))       mlm_chnk_size       <- 100
if (!exists("sweep_pct"))           sweep_pct           <- c(25, 5, 1, 0.2)
sweep_z <- round(qnorm(1 - sweep_pct / 100), 3)

## SETUP      -- 2x2 CELL GRID PARAMETERS --------------------------------------

# Allocation mechanism x variance structure: the paper's organising framework.
# variance_regime : homoscedastic   = constant variance across ages (1.0/1.0/1.0)
#                   heteroscedastic = rising   variance across ages (1.0/1.5/2.0)
# allocation_rule : random   = group assigned independent of baseline
#                   targeted = group assigned on baseline value (truncation)
# simulate        : TRUE -> fresh Monte Carlo run this build
# worked_example  : consistent worked example threaded throughout Intro /
#                   Methods / Discussion; B = ad-ban; D = letter

cell_grid <- data.frame(
  cell            = c("A",       "B",         "C",         "D"),
  variance_regime = c("homo",    "hetero",    "homo",      "hetero"),
  allocation_rule = c("random",  "random",    "targeted",  "targeted"),
  label           = c("Positive control / calibration check",
                      "Pop-wide measure (sugar levy / ad-ban)",
                      "Extreme Beggs et al.",
                      "Targeted baseline (\"letter\")"),
  simulate        = c(TRUE,      TRUE,        TRUE,        TRUE),
  worked_example  = c(NA,        "ad-ban",    NA,          "letter"),
  stringsAsFactors = FALSE)

# cache_version : manual cache-invalidation lever, folded into every cache_tag() 
# call (pools, MLM, comparison, sweep, diagnostic). The string's content carries 
# no information the script reads - only WHETHER it changed since the cache on
# disk was written matters. Bump it any changes are made to the population 
# generation, allocation, or replication-runner logic (including file name 
# construction logic) and every existing cache is invalidated, rather than
# deleting Data/ files by hand or hunting down which ones are stale. Bumping 
# does NOT delete old cache files - they are orphaned & can be cleaned up 
# separately if disk space matters.

cache_version <- "v4"

## SETUP      -- SCENARIO & SIMULATION SIZES -----------------------------------

# Ages at which weight is observed (years); N observations per cell; K (set under 
# RUN-CONTROL FLAGS above) is the shared replication count for every # Monte 
# Carlo pass in the script; ide_types crossed against every cell

ages        <- c(1, 2.5, 4)
Nsmp        <- 500
ide_types   <- c("no_ide", "true_ide")

# mean_mult: scales the POPULATION follow-up mean level (ages 2.5 & 4) for the 
# treatment arm, so the intervention age4 mean lands at 14.17 kg (1.5 kg/yr from 
# the 9.67 kg baseline over 3 yrs), against the control's ~2.07 kg/yr; written
# as 14.17/15.88 to keep the target legible. Baked into the treatment-arm 
# generating distribution, not applied afterwards.

# NOTE ON UNITS: sd_homo/sd_hetero/sd_mult_true_ide are kept in SD (not variance) 
# throughout script because mvrnorm() & covariance-matrix builder r2cov() -  
# Sigma = D %*% R %*% D, D = diag(sd) - both need SDs directly - expressed as 
# variances would require a sqrt() at every call with no benefit. Variance is 
# still the quantity the manuscript  any reporting discusses; SD is the internal
# computational unit for building Sigma.

# sd_mult_true_ide: the raw (marginal) SD ratio (vs control) reached at age4
# under true_ide, applied under a fixed correlation matrix - PRIMARY definition,
# build_treatment_distribution() - so the intervention weakens tracking (the
# absolute baseline-follow-up covariance) while leaving Corr(age1,age4)
# unchanged; the age2.5 ratio is interpolated proportional to elapsed time.
# The SAME parameter is reused, under a different interpretation (a
# CONDITIONAL-SD ratio given baseline), by the alternative Appendix-only
# definition in build_treat_dist_condshrink() - see PART TWO.

mean_mult        <- 14.17 / 15.88
sd_mult_true_ide <- 1.75  / 2.0

# All charts are 27 cm wide & 17 cm tall; multi-row calls pass plot_h explicitly

plot_w <- 10.63  ## 27 cm
plot_h <-  6.69  ## 17 cm

## SETUP      -- PROJECT ENVIRONMENT, PATHS, & DIRECTORIES ---------------------

# The working directory & the project '.Rproj' RStudio file name should be equal

check_project_name_consistency <- function() {
  wd_path     <- getwd()
  wd_name     <- basename(wd_path)
  rproj_files <- list.files(path = wd_path, pattern = "\\.Rproj$", 
                            full.names = FALSE)
  if (length(rproj_files) == 0) {
    message("FAIL  No .Rproj file found in: ", wd_path)
    return(invisible(FALSE)) }
  if (length(rproj_files) > 1) {
    message("WARN  Multiple .Rproj files found: ", 
            paste(rproj_files, collapse = ", "))
    message("      Expected exactly one.") }
  rproj_name <- tools::file_path_sans_ext(rproj_files[1])
  match      <- identical(wd_name, rproj_name)
  if (!match) {
    cat("Working directory name : ", wd_name,   "\n", sep = "")
    cat(".Rproj file name       : ", rproj_name, "\n", sep = "")
    cat("-------------------------------------------\n") }
  if (match) { message("Working directory is correct") } else {
    message("Mismatch between project name and working directory detected!") }
  return(invisible(match)) }
check_project_name_consistency()

# Set up storage directories, rooted at the working directory

setup_storage <- function() {
  wd_path       <- getwd()
  subdirs       <- c("Data", "Plots")
  ProjDir       <- wd_path
  SubDir        <- wd_path
  for (d in subdirs) { dir.create(file.path(SubDir, d), 
                                  showWarnings = FALSE, recursive = TRUE) }
  message("Storage rooted at: ", SubDir)
  paths         <- setNames(lapply(subdirs, function(d) file.path(SubDir, d)), 
                            subdirs)
  paths$ProjDir <- ProjDir
  return(invisible(paths)) }
storage_paths   <- setup_storage()
Data            <- storage_paths$Data
Plots           <- storage_paths$Plots
ProjDir         <- storage_paths$ProjDir

# All cache & results for this run live under Data/K_<K>/, not directly under
# Data/ - reassigning Data itself here (rather than editing every individual
# file.path(Data, ...) call site throughout the script) means every later
# stage - pools, verify/MLM results, appendix, checkpoints
# (run_reps_chunked()'s checkpoint_dir references Data via global
# scope), & the results xlsx - automatically lands under the K-specific
# subfolder with no other code changes needed. Lets you manually delete an
# entire K_<K> folder to drop all cache for that K value without touching any
# other K's cache, since nothing is shared across K subfolders.

Data_root       <- Data   # top-level /Data
Data            <- file.path(Data_root, sprintf("K_%d", K))
dir.create(Data, showWarnings = FALSE, recursive = TRUE)
message("Cache & results for this run rooted at: ", Data)

# Sweep_Data: base directory for the threshold-severity sweep (PART FOUR, both
# legs) & PART FIVE's diagnostic. Its cache follows the same K-subfolder rule 
# as everything else, Data/K_<K>/; sweep_K (PART FOUR) is folded into its cache 
# tags, so a sweep run under a different sweep_K_max sits beside it without 
# clashing.

Sweep_Data      <- Data

## SETUP      -- POPULATION MOMENTS (CDC/LMS) ----------------------------------

# Target means at ages 1 / 2.5 / 4; fixed population parameters, drawn from 
# CDC/LMS percentile data files (girls' weight-for-age)

mu_age1  <- 9.67
mu_age25 <- 12.78
mu_age4  <- 15.88

# Target standard deviations at ages 1 / 2.5 / 4 
# (kept as SD, not variance, for same reason noted under sd_mult_true_ide above: 
# r2cov() needs SDs to build Sigma)
# homoscedastic   regime: constant SD, set to the age-1 value throughout - 
#                         unrealistic but necessary to calibrate the evaluation 
#                         machinery
# heteroscedastic regime: rising SD, matching the original CDC/LMS-derived 
#                         spread - more naturally reflecting reality for growth 
#                         measures

sd_homo   <- c(1.0, 1.0, 1.0)
sd_hetero <- c(1.0, 1.5, 2.0)

# Lag autocorrelation: r = 0.9^delta_t between any two observed ages

rho_base <- 0.9

## #############################################################################
## PART TWO   -- FUNCTIONS LIBRARY
## FUNCTIONS  -- RANDOM SEED GENERATOR, LOGGING, & TIMING ----------------------

# Deterministic seed generator: polynomial rolling hash in modulus M = 2^31 - 1
# All arithmetic is performed in double precision to avoid integer overflow,
# with results truncated to integer only at the final step

seed_for <- function(...) {
  key    <- paste(..., collapse = "|")
  bytes  <- as.double(utf8ToInt(key))
  M      <- 2147483647.0   # 2^31 - 1 as double
  h      <- 0.0
  for (b in bytes) h <- (h * 131.0 + b) %% M
  h      <- as.integer(h)
  if (is.na(h) || h == 0L) h <- 1L
  return(h) }

# Parameter-hash cache tag: same rolling-hash idiom as seed_for(), but takes the
# actual PARAMETER VALUES feeding a cached object (not just labels) & returns a 
# short hex tag for the cache filename. Any change to a parameter passed in will 
# change the tag, so a stale cache is never loaded - a changed tag simply misses 
# the old file & regenerates under a new name

cache_tag <- function(...) {
  key    <- paste(..., collapse = "|")
  bytes  <- as.double(utf8ToInt(key))
  M      <- 2147483647.0   # 2^31 - 1 as double
  h      <- 0.0
  for (b in bytes) h <- (h * 131.0 + b) %% M
  h      <- as.integer(h)
  if (is.na(h) || h == 0L) h <- 1L
  return(sprintf("%08x", h)) }

# Multiple time stamp & logging routines

timestamp <- function() return(format(Sys.time(), "%H:%M:%S"))

format_seconds <- function(sec) {
  if (sec < 60)   return(sprintf("%.1f s",   sec))
  if (sec < 3600) return(sprintf("%.1f min", sec / 60))
  return(sprintf("%.1f hr", sec / 3600)) }

log_done  <- function(t0, label) {
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  log_step(paste0(label, " done  [", format_seconds(elapsed), "]"))
  return(invisible(NULL)) }

log_step  <- function(msg, Rows = 0) {
  if (Rows > 0) rep(cat("\n"), Rows)
  cat("[", timestamp(), "] ", msg, "\n", sep = "")
  return(invisible(NULL)) }

## FUNCTIONS  -- COMMON POPULATION GENERATOR -----------------------------------

# Converts a vector of lower-triangle correlations (read column-wise, e.g.
# c(r12, r13, r23) for a 3x3 matrix) into a full symmetric correlation matrix
# with unit diagonal

lower2R <- function(rvec) {
  k               <- (1 + sqrt(1 + 8 * length(rvec))) / 2
  R               <- diag(k)
  R[lower.tri(R)] <- rvec
  R[upper.tri(R)] <- t(R)[upper.tri(R)]
  return(R) }

# Converts a correlation matrix R & a vector of standard deviations sd into a
# covariance matrix: S = D %*% R %*% D, where D = diag(sd)

r2cov <- function(sd, R) {
  D   <- diag(sd)
  S   <- D %*% R %*% D
  return(S) }

# Draws the shared underlying population for one replication: Nobs children's
# weight trajectories at ages 1 / 2.5 / 4, multivariate normal, lag correlation
# r = rho_base^delta_age & population SD vector (kept as SD, not variance, for
# the reason noted under sd_mult_true_ide in PART ONE) chosen by the
# caller (sd_homo or sd_hetero). This is the OFAT anchor; the same population is
# reused across all four cells within a replication, so cross-cell differences 
# are attributable to the toggled factor alone - empirical = FALSE preserves 
# finite-sample sampling variation at N = 500, which the power/Type 1 figures 
# depend on. 

# Sigma_override: normally generate_population() builds Sigma from sd_vec &
# rho_base alone, via the fixed lag-decay recipe r2cov(sd_vec, R) where R's
# off-diagonal entries are all rho_base^delta_age - a single SD vector plus a
# single decay rate. Call sites pass the CONTROL sd_vec (sd_homo/sd_hetero) as
# the sd_vec argument while separately supplying the treatment arm's actual
# Sigma via Sigma_override, so an override is always needed mechanically, even
# though build_treatment_distribution()'s PRIMARY true_ide branch (rescale
# every occasion's raw SD under a fixed R) technically COULD be re-expressed
# via that same standard recipe with a different SD vector - it's supplied as
# an override rather than restructuring every call site's sd_vec argument.
# The ALTERNATIVE (Appendix-only) true_ide branch, build_treatment_
# distribution_condshrink(), genuinely CANNOT be expressed that way: it
# shrinks the conditional (post-baseline) variance directly, via partitioned-
# matrix algebra (Sigma11/Sigma21/Sigma22 blocks), then converts back to an
# unconditional covariance matrix - the result is a matrix shape no single
# SD vector under the fixed-decay R recipe can produce. Sigma_override is the
# escape hatch that lets generate_population() accept that pre-built matrix
# directly (see is.null(Sigma_override) below) rather than trying, and
# failing, to re-derive it from a plain SD vector & decay rate.

generate_population <- function(N, mu_vec, sd_vec, rho_base, seed, 
                                Sigma_override = NULL) {
  set.seed(seed)
  if (is.null(Sigma_override)) {
    delta_12  <- ages[2] - ages[1]
    delta_23  <- ages[3] - ages[2]
    delta_13  <- ages[3] - ages[1]
    rho_12    <- rho_base ^ delta_12
    rho_23    <- rho_base ^ delta_23
    rho_13    <- rho_base ^ delta_13
    Sigma     <- r2cov(sd = sd_vec, R = lower2R(c(rho_12, rho_13, rho_23))) 
  } else { Sigma <- Sigma_override }
  wt_mat      <- mvrnorm(N, mu_vec, Sigma, empirical = FALSE)
  pop         <- data.frame(
    id        = seq_len(N),
    wt_age1   = wt_mat[, 1],
    wt_age25  = wt_mat[, 2],
    wt_age4   = wt_mat[, 3],
    stringsAsFactors = FALSE)
  return(pop) }

# Compares a drawn population's sample moments (means, SDs, lag correlations)
# against the CDC/LMS population targets; returns a list of two data frames
# (moments, correlations). Used only in PART THREE's "VERIFY POPULATION
# GENERATOR" step, where its output is print()'d straight to the console for
# one homo & one hetero check population - no file is written & no MLwiN call
# is involved, so this check runs (and can be read) well before PART FOUR's
# MLM pass starts

verify_population <- function(pop, mu_vec, sd_vec, rho_base) {
  delta_12   <- ages[2] - ages[1]
  delta_23   <- ages[3] - ages[2]
  delta_13   <- ages[3] - ages[1]
  wt_cols    <- c("wt_age1", "wt_age25", "wt_age4")
  mean_hat   <- sapply(pop[, wt_cols], mean)
  sd_hat     <- sapply(pop[, wt_cols], sd)
  cor_mat    <- cor(pop[, wt_cols])
  out        <- data.frame(
    age       = ages,
    mean_hat  = round(as.numeric(mean_hat), 2),
    mean_tgt  = round(mu_vec, 2),
    sd_hat    = round(as.numeric(sd_hat), 2),
    sd_tgt    = round(sd_vec, 2),
    stringsAsFactors = FALSE)
  rho_check  <- data.frame(
    pair      = c("age1-age2.5", "age2.5-age4", "age1-age4"),
    rho_hat   = round(c(cor_mat[1, 2], cor_mat[2, 3], cor_mat[1, 3]), 2),
    rho_tgt   = round(c(rho_base ^ delta_12, rho_base ^ delta_23, 
                        rho_base ^ delta_13), 2),
    stringsAsFactors = FALSE)
  return(list(moments = out, correlations = rho_check)) }

## FUNCTIONS  -- TREATMENT-ARM DISTRIBUTION (IDE BAKED INTO GENERATION) --------

# Computes the treatment arm's generating mu_vec & Sigma for a given ide_type,
# baking the intervention effect directly into the distribution the treatment
# population is DRAWN from - not applied afterwards as a reconstruction.
# age1 is left at the control mean (the intervention has no baseline effect);
# ages 2.5 & 4 have their mean scaled by mean_mult.

# ide_type = "no_ide"  : mean shift only; Sigma unchanged from the control
#                        regime, so the raw (marginal) SD of every occasion,
#                        & hence baseline's predictive relationship with
#                        follow-up, is identical to control
# ide_type = "true_ide": as no_ide, but PRIMARY DEFINITION (used everywhere in
#                        the main analysis): the intervention weakens
#                        TRACKING itself - the biological/behavioural process
#                        by which early weight predicts later weight. Modelled
#                        as a uniform rescaling of the raw (marginal) SD at
#                        every follow-up occasion, reaching sd_mult_true_ide by
#                        age4 (age2.5 ratio interpolated proportional to
#                        elapsed time), under the SAME correlation matrix as
#                        control - since Cov = SD1*SD2*rho, this necessarily
#                        shrinks Cov(baseline, follow-up) right along with the
#                        marginal variances (Corr(baseline,age4) is therefore
#                        UNCHANGED from control - a pure "rescale everything,
#                        same shape" transform), which is what makes baseline
#                        a less powerful predictor of follow-up in absolute
#                        (covariance) terms despite the correlation itself
#                        being preserved. See build_treatment_distribution_
#                        condshrink() below for the alternative ("baseline's
#                        predictive influence protected") definition, kept
#                        only for the Appendix robustness check.

# Generating the treatment population this way means truncating on baseline
# afterwards - build_targeted_pool() - is applied identically to treatment &
# control pools, with no asymmetry between them

build_treatment_distribution <- function(mu_vec, sd_vec, rho_base, mean_mult, 
                                         sd_mult_true_ide, ide_type) {
  R_full       <- lower2R(c(
    rho_base ^ (ages[2] - ages[1]), 
    rho_base ^ (ages[3] - ages[1]), 
    rho_base ^ (ages[3] - ages[2])))
  Sigma_full   <- r2cov(sd = sd_vec, R = R_full)
  mu_int       <- c(mu_vec[1], mu_vec[2:3] * mean_mult)
  if (ide_type == "no_ide") { 
    return(list(mu_vec = mu_int, Sigma = Sigma_full)) }
  delta_25     <- ages[2] - ages[1]
  delta_4      <- ages[3] - ages[1]
  ratio_25     <- 1 - (1 - sd_mult_true_ide) * (delta_25 / delta_4)
  sd_int       <- sd_vec * c(1, ratio_25, sd_mult_true_ide)
  Sigma_int    <- r2cov(sd = sd_int, R = R_full)
  return(list(mu_vec = mu_int, Sigma = Sigma_int)) }

# ALTERNATIVE ("tracking-preserving") true_ide definition - APPENDIX ONLY.
# Rather than rescaling the raw marginal SD (& thus the raw baseline-follow-up
# covariance) under a fixed correlation matrix, this instead partitions Sigma
# into the component EXPLAINED by baseline (Sigma21 %*% t(Sigma21) / Sigma11,
# i.e., the regression-on-baseline part) & the CONDITIONAL/residual component
# left over given baseline, then shrinks ONLY the conditional component by
# sd_mult_true_ide, leaving the raw baseline-follow-up covariance (Sigma21)
# untouched. This is the definition the main script used before the switch
# to build_treatment_distribution() above; it represents an intervention
# that adds consistency around each child's OWN baseline-predicted
# trajectory without touching the tracking relationship itself (Corr(baseline,
# age4) is NOT preserved under this version - it slightly INCREASES, since
# the numerator covariance is fixed while total follow-up variance shrinks).
# Used only to build a supplementary Appendix Table 2 (Oldham + MLM, cells
# A-D, Type 1 error & power) as a robustness/sensitivity check against the
# primary tracking-weakening assumption - NOT part of the main-text pipeline.

build_treat_dist_condshrink <- function(mu_vec, sd_vec, rho_base, 
                                        mean_mult, sd_mult_true_ide, ide_type) {
  Sigma_full   <- r2cov(sd = sd_vec, R = lower2R(c(
    rho_base ^ (ages[2] - ages[1]), 
    rho_base ^ (ages[3] - ages[1]), 
    rho_base ^ (ages[3] - ages[2]))))
  mu_int       <- c(mu_vec[1], mu_vec[2:3] * mean_mult)
  if (ide_type == "no_ide") { return(list(mu_vec = mu_int, Sigma = Sigma_full)) }
  Sigma11      <- Sigma_full[1, 1]
  Sigma21      <- Sigma_full[2:3, 1]
  Sigma22      <- Sigma_full[2:3, 2:3]
  Sigma_cond   <- Sigma22 - outer(Sigma21, Sigma21) / Sigma11
  delta_25     <- ages[2] - ages[1]
  delta_4      <- ages[3] - ages[1]
  ratio_25     <- 1 - (1 - sd_mult_true_ide) * (delta_25 / delta_4)
  sd_cond      <- sqrt(diag(Sigma_cond))
  R_cond       <- Sigma_cond / outer(sd_cond, sd_cond)
  sd_cond_neg  <- sd_cond * c(ratio_25, sd_mult_true_ide)
  Sigma_cond_neg <- outer(sd_cond_neg, sd_cond_neg) * R_cond
  Sigma_int    <- Sigma_full
  Sigma_int[2:3, 2:3] <- Sigma_cond_neg + outer(Sigma21, Sigma21) / Sigma11
  return(list(mu_vec = mu_int, Sigma = Sigma_int)) }

## FUNCTIONS  -- TARGETED-ALLOCATION POOLS (CELLS C & D) -----------------------

# Targeted cells (C, D) do not threshold a 500-row sample; in reality the
# treatment & control arms are different children, so each arm is drawn
# independently from the at-risk tail of its own large reference population.
# Both arms are thresholded identically: the control arm must be drawn from the 
# SAME restricted / truncated population as the treatment arm, since an
# exchangeable comparator has to share the treatment arm's baseline variance 
# structure - an unthresholded control would have a normal baseline variance
# while the treatment arm's is truncation-compressed, which breaks Oldham's 
# calibration regardless of any IDE effect.

# A "pool" is the SPECIFICATION of one arm's truncated distribution: the
# untruncated mean vector & covariance matrix under the given variance regime,
# plus the baseline threshold (the at-risk cut, e.g., z > 1.645 against the
# CDC/LMS reference). Every replication draws afresh from this exact
# distribution rather than from one finite pre-built population: a fixed
# finite population carries its own sampling error, shared by all K
# replications, which does not shrink as K grows & would move power by
# several percentage points at the more severe thresholds

build_targeted_pool <- function(mu_vec, sd_vec, rho_base, ref_mean, ref_sd,
                                Sigma_override = NULL, z_cut = 1.645) {
  if (is.null(Sigma_override)) {
    Sigma      <- r2cov(sd = sd_vec, R = lower2R(c(
      rho_base ^ (ages[2] - ages[1]), 
      rho_base ^ (ages[3] - ages[1]), 
      rho_base ^ (ages[3] - ages[2]))))
  } else { Sigma <- Sigma_override }
  threshold    <- ref_mean[1] + z_cut * ref_sd[1]
  return(list(mu_vec = mu_vec, Sigma = Sigma, threshold = threshold)) }

# Draws exactly n_draw (= Nsmp = 500) children for one replication from the
# exact baseline-truncated distribution in pool (build_targeted_pool() above):
# baseline weight from the normal distribution truncated below at the
# threshold (inverse CDF on the upper tail, stable even at the 0.2% cut), then
# weights at ages 2.5 & 4 from their normal distribution conditional on that
# baseline. Truncating on baseline leaves the conditional distribution of the
# later ages unchanged, so this reproduces exactly a multivariate normal
# population thresholded on baseline, with no finite pool behind it

sample_from_pool <- function(pool, n_draw, seed) { set.seed(seed)
  mu             <- pool$mu_vec
  S              <- pool$Sigma
  sd1            <- sqrt(S[1, 1])
  tail_p         <- pnorm((pool$threshold - mu[1]) / sd1, lower.tail = FALSE)
  wt_age1        <- mu[1] + sd1 * qnorm(runif(n_draw) * tail_p, 
                                        lower.tail = FALSE)
  beta           <- S[2:3, 1] / S[1, 1]
  S_cond         <- S[2:3, 2:3] - outer(S[2:3, 1], S[2:3, 1]) / S[1, 1]
  later          <- outer(wt_age1 - mu[1], beta) + mvrnorm(n_draw, mu[2:3], 
                                                           S_cond)
  out            <- data.frame(
    id           = seq_len(n_draw),
    wt_age1      = wt_age1,
    wt_age25     = later[, 1],
    wt_age4      = later[, 2],
    stringsAsFactors = FALSE)
  return(out) }

## FUNCTIONS  -- OLDHAM'S METHOD -----------------------------------------------

# Oldham's correlation for one group: correlates the difference between baseline 
# & follow-up - age1-age4 - with their average - (age1+age4)/2 - which is 
# uncorrelated by construction when variances are equal. Returns the correlation, 
# its cor.test() against zero, & the Fisher z-transform

oldham_one_group <- function(wt_age1, wt_age4) {
  diff_val  <- wt_age1 - wt_age4
  ave_val   <- (wt_age1 + wt_age4) / 2
  test_obj  <- cor.test(diff_val, ave_val, alternative = "two.sided", 
                        method = "pearson")
  r_val     <- as.numeric(test_obj$estimate)
  z_val     <- fisherz(r_val)
  return(list(r = r_val, p = test_obj$p.value, z = z_val, n = length(wt_age1))) }

# Run Oldham's method on one allocated & IDE-injected population, returning both 
# contrasts used throughout the paper:
#   intervention_only : group1's Oldham correlation tested against zero; the 
#                       naive single-group test fails under heteroscedasticity 
#                       (Cells B-D) even with no IDE
#   vs_comparator     : group1 vs group0's Oldham correlations, contrasted via 
#                       Fisher's z; a valid counter-factual for Cell B (random 
#                       allocation, where the contrast restores calibration) but 
#                       only with a like-for-like comparator

run_oldham <- function(pop) {
  g1     <- pop[pop$group == 1, ]
  g0     <- pop[pop$group == 0, ]
  old_1  <- oldham_one_group(g1$wt_age1, g1$wt_age4)
  old_0  <- oldham_one_group(g0$wt_age1, g0$wt_age4)
  z_test <- cocor.indep.groups(old_0$r, old_1$r, old_0$n, 
                               old_1$n, alternative = "two.sided", 
                               test = "fisher1925", return.htest = TRUE)
  intervention_only <- list(r = old_1$r, p = old_1$p)
  vs_comparator     <- list(r_comparator = old_0$r, 
                            r_intervention = old_1$r, 
                            p = z_test$fisher1925$p.value)
  return(list(intervention_only = intervention_only, 
              vs_comparator = vs_comparator)) }

## FUNCTIONS  -- ONE REPLICATION, ALL FOUR CELLS (OLDHAM ROUTE) ----------------

# One replication of cell A/B/C/D under BOTH no_ide & true_ide, returning each
# cell's intervention-only & (where applicable) vs-comparator p-values for both.
# no_ide p-values feed Type 1 error; true_ide p-values feed power.

# pool_set: named list with treat_C, control_C, treat_D, control_D, plus the
# true_ide treatment pools treat_C_neg & treat_D_neg - the truncated pools to
# sample from for the targeted cells. Defaults to the primary 5% threshold pool.
# dist_A_neg_arg / dist_B_neg_arg: the true_ide distributions for cells A & B,
# defaulting to the primary (tracking-weakening) globals dist_A_true_ide /
# dist_B_true_ide - overridable so the SAME replication function can be reused
# unchanged for the Appendix run against the alternative (tracking-preserving)
# distributions, without touching the primary code path. no_ide is identical
# under either definition (see build_treatment_distribution()), so dist_A_no_ide
# / dist_B_no_ide are not parameterised - always the primary globals.
# Relies on dist_A_no_ide / dist_B_no_ide & the six pool_* objects existing in
# the calling environment before use - both built in PART FOUR, ahead of the
# first call to this function.

run_one_replication <- function(rep_i, pool_set = list(
  treat_C     = pool_treat_C_no_ide,
  treat_C_neg = pool_treat_C_true_ide,
  control_C   = pool_control_C,
  treat_D     = pool_treat_D_no_ide,
  treat_D_neg = pool_treat_D_true_ide,
  control_D   = pool_control_D),
  dist_A_neg_arg = dist_A_true_ide, dist_B_neg_arg = dist_B_true_ide) {
  pop_A        <- generate_population(N = Nsmp, mu_vec = dist_A_no_ide$mu_vec, 
                                      sd_vec = sd_homo, rho_base = rho_base, 
                                      seed = seed_for("verify", rep_i, 
                                                      "pop", "A"), 
                                      Sigma_override = dist_A_no_ide$Sigma)
  pop_A_neg    <- generate_population(N = Nsmp, mu_vec = dist_A_neg_arg$mu_vec, 
                                      sd_vec = sd_homo, rho_base = rho_base, 
                                      seed = seed_for("verify", rep_i, 
                                                      "pop", "A", "neg"), 
                                      Sigma_override = dist_A_neg_arg$Sigma)
  pop_A$group  <- 1; pop_A_neg$group <- 1
  draw_treat_B <- generate_population(N = Nsmp, mu_vec = dist_B_no_ide$mu_vec, 
                                      sd_vec = sd_hetero, rho_base = rho_base, 
                                      seed = seed_for("verify", rep_i, "pop", 
                                                      "B", "treat"), 
                                      Sigma_override = dist_B_no_ide$Sigma)
  draw_treat_B_neg <- generate_population(N = Nsmp, 
                                          mu_vec = dist_B_neg_arg$mu_vec, 
                                          sd_vec = sd_hetero, 
                                          rho_base = rho_base, 
                                          seed = seed_for("verify", rep_i, 
                                                          "pop", "B", "treat", 
                                                          "neg"),
                                          Sigma_override = dist_B_neg_arg$Sigma)
  draw_control_B   <- generate_population(N = Nsmp, mu_vec = mu_pop, 
                                          sd_vec = sd_hetero, 
                                          rho_base = rho_base, 
                                          seed = seed_for("verify", rep_i, 
                                                          "pop", "B", 
                                                          "control"))
  draw_treat_B$group     <- 1; draw_control_B$group <- 0
  draw_treat_B_neg$group <- 1
  draw_treat_C           <- sample_from_pool(pool_set$treat_C,     
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", rep_i, 
                                                             "sample", "treat",
                                                             "C"))
  draw_treat_C_neg       <- sample_from_pool(pool_set$treat_C_neg, 
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", rep_i, 
                                                             "sample", "treat",
                                                             "C", "neg"))
  draw_control_C         <- sample_from_pool(pool_set$control_C,   
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", rep_i, 
                                                             "sample", 
                                                             "control", "C"))
  draw_treat_D           <- sample_from_pool(pool_set$treat_D,     
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", 
                                                             rep_i, "sample", 
                                                             "treat",   "D"))
  draw_treat_D_neg       <- sample_from_pool(pool_set$treat_D_neg, 
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", rep_i, 
                                                             "sample", "treat",
                                                             "D", "neg"))
  draw_control_D         <- sample_from_pool(pool_set$control_D,   
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", rep_i, 
                                                             "sample", 
                                                             "control", "D"))
  draw_treat_C$group     <- 1; draw_control_C$group <- 0
  draw_treat_C_neg$group <- 1
  draw_treat_D$group     <- 1; draw_control_D$group <- 0
  draw_treat_D_neg$group <- 1
  cA      <- pop_A
  cA_neg  <- pop_A_neg
  cB      <- rbind(draw_control_B, draw_treat_B)
  cB_neg  <- rbind(draw_control_B, draw_treat_B_neg)
  cC      <- rbind(draw_control_C, draw_treat_C)
  cC_neg  <- rbind(draw_control_C, draw_treat_C_neg)
  cD      <- rbind(draw_control_D, draw_treat_D)
  cD_neg  <- rbind(draw_control_D, draw_treat_D_neg)
  oA      <- oldham_one_group(cA$wt_age1, cA$wt_age4)
  oA_neg  <- oldham_one_group(cA_neg$wt_age1, cA_neg$wt_age4)
  oB      <- run_oldham(cB)
  oB_neg  <- run_oldham(cB_neg)
  oC      <- run_oldham(cC)
  oC_neg  <- run_oldham(cC_neg)
  oD      <- run_oldham(cD)
  oD_neg  <- run_oldham(cD_neg)
  return(data.frame(
    rep                   = rep_i,
    p_A_int_only          = oA$p,
    p_A_int_only_neg      = oA_neg$p,
    p_B_int_only          = oB$intervention_only$p,
    p_B_vs_comparator     = oB$vs_comparator$p,
    p_B_int_only_neg      = oB_neg$intervention_only$p,
    p_B_vs_comparator_neg = oB_neg$vs_comparator$p,
    p_C_int_only          = oC$intervention_only$p,
    p_C_vs_comparator     = oC$vs_comparator$p,
    p_C_int_only_neg      = oC_neg$intervention_only$p,
    p_C_vs_comparator_neg = oC_neg$vs_comparator$p,
    p_D_int_only          = oD$intervention_only$p,
    p_D_vs_comparator     = oD$vs_comparator$p,
    p_D_int_only_neg      = oD_neg$intervention_only$p,
    p_D_vs_comparator_neg = oD_neg$vs_comparator$p) ) }

# Power-only variant of run_one_replication(): computes ONLY the true_ide
# (power-side) p-values, skipping the no_ide draws & tests entirely. Used
# ONLY by the Appendix run (PART FOUR) - the Appendix's Type 1 error figures
# are mathematically guaranteed to match the primary run's (no_ide is
# unaffected by which true_ide builder is used, see build_treatment_
# distribution()), so they are reused directly from verify_res rather than
# recomputed, roughly halving the Appendix Oldham pass's cost. Control-arm
# draws are STILL needed (shared with the true_ide comparator test) & so are
# NOT skipped - only the no_ide TREATMENT draws & the no_ide tests are
# dropped. Seed keys for every retained draw are identical to
# run_one_replication()'s, so results are byte-for-byte the same as the
# true_ide half of a full run_one_replication() call - this is purely an
# efficiency change, not a different computation.

run_one_rep_power_only <- function(rep_i, pool_set = list(
  treat_C_neg = pool_treat_C_true_ide,
  control_C   = pool_control_C,
  treat_D_neg = pool_treat_D_true_ide,
  control_D   = pool_control_D),
  dist_A_neg_arg = dist_A_true_ide, dist_B_neg_arg = dist_B_true_ide) {
  pop_A_neg        <- generate_population(N = Nsmp, 
                                          mu_vec = dist_A_neg_arg$mu_vec, 
                                          sd_vec = sd_homo,
                                          rho_base = rho_base, 
                                          seed = seed_for("verify", rep_i, 
                                                          "pop", "A", "neg"),
                                          Sigma_override = dist_A_neg_arg$Sigma)
  pop_A_neg$group  <- 1
  draw_treat_B_neg <- generate_population(N = Nsmp, 
                                          mu_vec = dist_B_neg_arg$mu_vec, 
                                          sd_vec = sd_hetero,
                                          rho_base = rho_base, 
                                          seed = seed_for("verify", rep_i, 
                                                          "pop", "B", "treat", 
                                                          "neg"),
                                          Sigma_override = dist_B_neg_arg$Sigma)
  draw_control_B   <- generate_population(N = Nsmp, mu_vec = mu_pop, 
                                          sd_vec = sd_hetero, 
                                          rho_base = rho_base, 
                                          seed = seed_for("verify", rep_i, 
                                                          "pop", "B", 
                                                          "control"))
  draw_treat_B_neg$group <- 1; draw_control_B$group <- 0
  draw_treat_C_neg       <- sample_from_pool(pool_set$treat_C_neg, 
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", rep_i, 
                                                             "sample", "treat",
                                                             "C", "neg"))
  draw_control_C         <- sample_from_pool(pool_set$control_C,   
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", rep_i, 
                                                             "sample", "control",
                                                             "C"))
  draw_treat_D_neg       <- sample_from_pool(pool_set$treat_D_neg, 
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", rep_i, 
                                                             "sample", "treat",
                                                             "D", "neg"))
  draw_control_D         <- sample_from_pool(pool_set$control_D,   
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", rep_i, 
                                                             "sample", 
                                                             "control", "D"))
  draw_treat_C_neg$group <- 1; draw_control_C$group <- 0
  draw_treat_D_neg$group <- 1; draw_control_D$group <- 0
  cA_neg  <- pop_A_neg
  cB_neg  <- rbind(draw_control_B, draw_treat_B_neg)
  cC_neg  <- rbind(draw_control_C, draw_treat_C_neg)
  cD_neg  <- rbind(draw_control_D, draw_treat_D_neg)
  oA_neg  <- oldham_one_group(cA_neg$wt_age1, cA_neg$wt_age4)
  oB_neg  <- run_oldham(cB_neg)
  oC_neg  <- run_oldham(cC_neg)
  oD_neg  <- run_oldham(cD_neg)
  return(data.frame(
    rep                    = rep_i,
    p_A_int_only_neg       = oA_neg$p,
    p_B_int_only_neg       = oB_neg$intervention_only$p,
    p_B_vs_comparator_neg  = oB_neg$vs_comparator$p,
    p_C_int_only_neg       = oC_neg$intervention_only$p,
    p_C_vs_comparator_neg  = oC_neg$vs_comparator$p,
    p_D_int_only_neg       = oD_neg$intervention_only$p,
    p_D_vs_comparator_neg  = oD_neg$vs_comparator$p) ) }

# Version for the threshold sweep: only for C & D (A & B do not use threshold-
# specific pools so their results are invariant across the sweep & need not be
# recomputed for every threshold). pool_set is looked up from the global
# sweep_pools list by z_tag key

# Draws BOTH no_ide (feeds Type 1 error) & true_ide (feeds power) treatment arms
# against the same shared control draw, mirroring run_one_replication()'s
# no_ide/true_ide pairing (PART TWO) - so the sweep can report power at every
# threshold, not just Type 1 error.

run_one_replication_CD <- function(rep_i, z_tag) {
  pool_set           <- sweep_pools[[z_tag]]$pool_set
  draw_treat_C       <- sample_from_pool(pool_set$treat_C,     
                                         n_draw = Nsmp, 
                                         seed = seed_for("sweep", rep_i, 
                                                         "sample", "treat",
                                                         "C"))
  draw_treat_C_neg   <- sample_from_pool(pool_set$treat_C_neg, 
                                         n_draw = Nsmp, 
                                         seed = seed_for("sweep", rep_i, 
                                                         "sample", "treat",
                                                         "C", "neg"))
  draw_control_C     <- sample_from_pool(pool_set$control_C,   
                                         n_draw = Nsmp, 
                                         seed = seed_for("sweep", rep_i, 
                                                         "sample", "control",
                                                         "C"))
  draw_treat_D       <- sample_from_pool(pool_set$treat_D,     
                                         n_draw = Nsmp, 
                                         seed = seed_for("sweep", rep_i, 
                                                         "sample", "treat",
                                                         "D"))
  draw_treat_D_neg   <- sample_from_pool(pool_set$treat_D_neg, 
                                         n_draw = Nsmp, 
                                         seed = seed_for("sweep", rep_i, 
                                                         "sample", "treat",
                                                         "D", "neg"))
  draw_control_D     <- sample_from_pool(pool_set$control_D,   
                                         n_draw = Nsmp, 
                                         seed = seed_for("sweep", rep_i, 
                                                         "sample", "control",
                                                         "D"))
  draw_treat_C$group <- 1; draw_treat_C_neg$group <- 1; draw_control_C$group <- 0
  draw_treat_D$group <- 1; draw_treat_D_neg$group <- 1; draw_control_D$group <- 0
  cC      <- rbind(draw_control_C, draw_treat_C)
  cC_neg  <- rbind(draw_control_C, draw_treat_C_neg)
  cD      <- rbind(draw_control_D, draw_treat_D)
  cD_neg  <- rbind(draw_control_D, draw_treat_D_neg)
  oC      <- run_oldham(cC)
  oC_neg  <- run_oldham(cC_neg)
  oD      <- run_oldham(cD)
  oD_neg  <- run_oldham(cD_neg)
  return(data.frame(
    rep                    = rep_i,
    p_C_int_only           = oC$intervention_only$p,
    p_C_vs_comparator      = oC$vs_comparator$p,
    p_C_int_only_neg       = oC_neg$intervention_only$p,
    p_C_vs_comparator_neg  = oC_neg$vs_comparator$p,
    p_D_int_only           = oD$intervention_only$p,
    p_D_vs_comparator      = oD$vs_comparator$p,
    p_D_int_only_neg       = oD_neg$intervention_only$p,
    p_D_vs_comparator_neg  = oD_neg$vs_comparator$p)) }

# MLM analogue of run_one_replication_CD(), used by the threshold-severity sweep
# so Oldham & MLM are compared on IDENTICAL draws (same seed keys, same pools)
# at every threshold. Defined here (rather than deferred to where MLM functions
# live) so it sits next to its Oldham counterpart; it calls mlm_contrasts()
# (PART TWO) which itself depends on lavaan_one_arm()/lavaan_two_arm(), both
# already defined above this point in the script. Also draws both no_ide &
# true_ide treatment arms, for the same reason as run_one_replication_CD() above.
# Returns a conv_* flag beside every p_* column, as run_one_replication_mlm_CD()
# does, so mlm_exclusions() can count non-converged fits per threshold.

run_one_replication_CD_mlm <- function(rep_i, z_tag) {
  pool_set           <- sweep_pools[[z_tag]]$pool_set
  draw_treat_C       <- sample_from_pool(pool_set$treat_C,     
                                         n_draw = Nsmp, 
                                         seed = seed_for("sweep", rep_i, 
                                                         "sample", "treat",
                                                         "C"))
  draw_treat_C_neg   <- sample_from_pool(pool_set$treat_C_neg, 
                                         n_draw = Nsmp, 
                                         seed = seed_for("sweep", rep_i, 
                                                         "sample", "treat",
                                                         "C", "neg"))
  draw_control_C     <- sample_from_pool(pool_set$control_C,   
                                         n_draw = Nsmp, 
                                         seed = seed_for("sweep", rep_i, 
                                                         "sample", "control",
                                                         "C"))
  draw_treat_D       <- sample_from_pool(pool_set$treat_D,     
                                         n_draw = Nsmp, 
                                         seed = seed_for("sweep", rep_i, 
                                                         "sample", "treat",
                                                         "D"))
  draw_treat_D_neg   <- sample_from_pool(pool_set$treat_D_neg, 
                                         n_draw = Nsmp, 
                                         seed = seed_for("sweep", rep_i, 
                                                         "sample", "treat",
                                                         "D", "neg"))
  draw_control_D     <- sample_from_pool(pool_set$control_D,   
                                         n_draw = Nsmp, 
                                         seed = seed_for("sweep", rep_i,
                                                         "sample", "control",
                                                         "D"))
  draw_treat_C$group <- 1; draw_treat_C_neg$group <- 1; draw_control_C$group <- 0
  draw_treat_D$group <- 1; draw_treat_D_neg$group <- 1; draw_control_D$group <- 0
  cC      <- rbind(draw_control_C, draw_treat_C)
  cC_neg  <- rbind(draw_control_C, draw_treat_C_neg)
  cD      <- rbind(draw_control_D, draw_treat_D)
  cD_neg  <- rbind(draw_control_D, draw_treat_D_neg)
  mC      <- mlm_contrasts(cC)
  mC_neg  <- mlm_contrasts(cC_neg)
  mD      <- mlm_contrasts(cD)
  mD_neg  <- mlm_contrasts(cD_neg)
  return(data.frame(
    rep                      = rep_i,
    p_C_int_only             = mC$intervention_only$p,
    conv_C_int_only          = mC$intervention_only$converged,
    p_C_vs_comparator        = mC$vs_comparator$p,
    conv_C_vs_comparator     = mC$vs_comparator$converged,
    p_C_int_only_neg         = mC_neg$intervention_only$p,
    conv_C_int_only_neg      = mC_neg$intervention_only$converged,
    p_C_vs_comparator_neg    = mC_neg$vs_comparator$p,
    conv_C_vs_comparator_neg = mC_neg$vs_comparator$converged,
    p_D_int_only             = mD$intervention_only$p,
    conv_D_int_only          = mD$intervention_only$converged,
    p_D_vs_comparator        = mD$vs_comparator$p,
    conv_D_vs_comparator     = mD$vs_comparator$converged,
    p_D_int_only_neg         = mD_neg$intervention_only$p,
    conv_D_int_only_neg      = mD_neg$intervention_only$converged,
    p_D_vs_comparator_neg    = mD_neg$vs_comparator$p,
    conv_D_vs_comparator_neg = mD_neg$vs_comparator$converged)) }

## FUNCTIONS  -- RUN REPLICATIONS IN CHECKPOINTED CHUNKS -----------------------

# Runs K replicates in fixed-size batches of chunk_size, logging progress after
# each batch (chunked progress reporting). Serial lapply() below parallel_K;
# parLapply() above it, since cluster start-up overhead exceeds the benefit for
# small K.

run_reps_chunked <- function(K, parallel_K,
                                     chunk_size     = 500,
                                     rep_fn         = run_one_replication,
                                     cl             = NULL,
                                     checkpoint_tag = NULL,
                                     force_fresh    = FALSE,
                                     data_dir       = Data, ...) {
  use_parallel   <- K > parallel_K
  extra_args     <- list(...)
  # CHECKPOINTING: when checkpoint_tag is supplied, each completed chunk is
  # written to its OWN small file under <data_dir>/checkpoints/<tag>/, rather
  # than one growing file rewritten IN FULL after every chunk. The one-big-
  # file design used to mean a cold-start run's total checkpoint I/O grew
  # with the SQUARE of the chunk count (chunk 1 rewrites 1 chunk's worth,
  # chunk 2 rewrites 2 chunks' worth, ..., chunk n rewrites all n) - for a
  # K = 1,000,000, chunk_size = 100 run that is 10,000 rewrites of an
  # ever-larger object, which turned into a multi-hour I/O bottleneck even
  # though the underlying replications were cheap. Writing only the new
  # chunk keeps each write's cost fixed at chunk_size, however far into the
  # run this is. It also means a corrupted/truncated file from a kill mid-
  # write only costs that one chunk, not the whole accumulated history, as
  # the old single-file design risked. Every replication's seed is derived
  # purely from rep_i (seed_for()), never from execution order or chunk
  # position, so resuming & dispatching only the missing rep_ids (identified
  # from each chunk file's own rep column, so no separate rep_ids index is
  # needed) produces identical results to an uninterrupted run - making it
  # safe to resume. Granularity is chunk_size: a crash mid-chunk still loses
  # that chunk's replications, since they are only written once the whole
  # chunk returns. The force_fresh flag deletes any existing checkpoint
  # before starting, since a deliberate force-rerun means start clean, not
  # silently resume from old progress. data_dir defaults to the global Data
  # (the current run's K-specific subfolder) but is overridable - the
  # threshold-severity sweep's call sites pass Sweep_Data instead, so its
  # checkpoints land alongside its other cache files (Sweep_Data, PART ONE),
  # which is now the same K-subfolder.
  checkpoint_dir <- if (!is.null(checkpoint_tag) ) 
    file.path(data_dir, "checkpoints", checkpoint_tag) else NULL
  if (force_fresh && !is.null(checkpoint_dir) && dir.exists(checkpoint_dir) ) {
    unlink(checkpoint_dir, recursive = TRUE)
    log_step("force_rerun set -- discarding existing checkpoint & starting fresh") }
  done_results       <- list()
  done_rep_ids       <- integer(0)
  if (!is.null(checkpoint_dir) && dir.exists(checkpoint_dir) ) {
    chunk_files      <- list.files(checkpoint_dir, 
                                   pattern = "\\.rds$", full.names = TRUE)
    if (length(chunk_files) > 0) {
      done_results   <- lapply(chunk_files, readRDS)
      done_rep_ids   <- unlist(lapply(done_results, function(d) d$rep), 
                               use.names = FALSE)
      log_step(sprintf("resuming from checkpoint: %d of %d replications already done", 
                       length(done_rep_ids), K) ) } }
  remaining_ids      <- setdiff(seq_len(K), done_rep_ids)
  if (length(remaining_ids) == 0) {
    log_step("all replications already complete in checkpoint -- nothing to dispatch")
    return(do.call(rbind, done_results) ) }
  remaining_chunk_bounds <- unique(c(seq(0, length(remaining_ids), 
                                         by = chunk_size), length(remaining_ids)))
  if (use_parallel) {
    # cl is always the shared script-wide cluster from get_script_cl(), already
    # built with libraries loaded; just re-export the global environment, since
    # closures like z_tag-specific pool sets can change between calls that
    # share the same persistent cluster
    clusterExport(cl, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
    n_chunks   <- length(remaining_chunk_bounds) - 1
    log_step(sprintf("dispatching %d replications in %d %s of %d", 
                     length(remaining_ids), n_chunks,
                     if (n_chunks == 1) "chunk" else "chunks", chunk_size)) }
  if (!is.null(checkpoint_dir) ) dir.create(checkpoint_dir, showWarnings = FALSE, 
                                            recursive = TRUE)
  bar_width    <- 40
  for (chunk_i in seq_len(length(remaining_chunk_bounds) - 1)) {
    rep_lo     <- remaining_chunk_bounds[chunk_i] + 1
    rep_hi     <- remaining_chunk_bounds[chunk_i + 1]
    rep_ids    <- remaining_ids[rep_lo:rep_hi]
    chunk_res  <- if (use_parallel) { parLapply(cl, rep_ids, rep_fn, ...)
    } else { lapply(rep_ids, rep_fn, ...) }
    chunk_df   <- do.call(rbind, chunk_res)
    done_results[[length(done_results) + 1]] <- chunk_df
    done_rep_ids <- c(done_rep_ids, rep_ids)
    if (!is.null(checkpoint_dir) ) {
      chunk_file <- file.path(checkpoint_dir, 
                              sprintf("chunk %010d-%010d.rds", 
                                      min(rep_ids), max(rep_ids)))
      saveRDS(chunk_df, chunk_file) }
    pct        <- length(done_rep_ids) / K
    n_filled   <- round(bar_width * pct)
    bar        <- paste0("[", strrep("=", n_filled), 
                         strrep(" ", bar_width - n_filled), "]")
    cat(sprintf("\r  %s %6.2f%%  %s / %s", bar, 100 * pct,
                format(as.integer(length(done_rep_ids)), big.mark = ","), 
                format(as.integer(K), big.mark = ",")))
    flush.console() }
  cat("\n")
  final_res <- do.call(rbind, done_results)
  if (!is.null(checkpoint_dir) && dir.exists(checkpoint_dir) ) 
    unlink(checkpoint_dir, recursive = TRUE)
  return(final_res) }

## FUNCTIONS  -- SCRIPT-WIDE CLUSTER (BUILT ONLY WHEN FIRST NEEDED) ------------

# One cluster for the whole script run, shared across every stage, but NOT built
# pre-emptively. get_script_cl() builds the cluster on its FIRST call that
# actually needs parallelism (from whichever stage is first to find no cache &
# need to dispatch replications) & returns the same handle on every later call
# - a fully-cached run builds no cluster & a run needing only one stage's worth
# of computation pays the start-up cost exactly once, not once per stage & not
# once regardless of need.

# needs_parallel (per-call, defaults to the top-level K > parallel_K test):
# call sites that pass their OWN parallel_K to run_reps_chunked (MLM
# pass, MLM sweep leg - both force parallel_K = 1 to always dispatch in
# parallel regardless of top-level K) must also pass that same effective
# comparison here, or this function would judge "is a cluster needed" purely
# from the top-level K/parallel_K, decide "no" & cache that NULL answer
# permanently (script_cl_built latches on the FIRST call) - leaving every
# later stage that actually does need one to build & tear down its own
# private cluster instead of sharing this one. That mismatch was the actual
# cause of "starting cluster" logging twice when top-level K <= parallel_K:
# each MLM-type stage silently fell back to run_reps_chunked's own
# owns_cluster path. Once script_cl_built is TRUE the decision is locked in
# for the rest of the run (matching a single shared cluster's lifetime) -
# it does NOT re-evaluate needs_parallel on later calls.

script_cl_built    <- FALSE
script_cl_handle   <- NULL

get_script_cl <- function(needs_parallel = (K > parallel_K)) {
  if (script_cl_built) return(script_cl_handle)
  script_cl_built <<- TRUE
  if (!needs_parallel) {
    script_cl_handle <<- NULL
    return(NULL) }
  n_cores        <- detectCores()
  log_step(sprintf("STARTING SCRIPT-WIDE CLUSTER (%d cores) ", n_cores), 1)
  cl_tmp         <- makeCluster(n_cores)
  clusterExport(cl_tmp, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
  clusterEvalQ(cl_tmp, { library(MASS); library(cocor); library(psych); 
    library(lavaan); return(NULL) })
  log_step("script-wide cluster now ready")
  script_cl_handle <<- cl_tmp
  return(script_cl_handle) }

## FUNCTIONS  -- REJECTION-RATE HELPERS (TYPE 1 ERROR / POWER) -----------------

# Compute the proportion of p-values at or below 0.05 across a vector of reps.
# na.rm = TRUE excludes any replication whose p-value is NA (in practice, only
# possible for MLM/lavaan replications where their tryCatch might catch a fit
# error or non-convergence); Oldham's closed-form correlation never produces
# NA. Without na.rm, a single NA silently blanks the entire Type 1 error
# figure to NA. The rate is therefore computed over CONVERGED reps only

reject_rate <- function(p_vec) return(100 * mean(p_vec <= 0.05, na.rm = TRUE))

# power_rate() is reject_rate() applied to the true_ide p-value vector, kept as a 
# distinct name (rather than reusing reject_rate() under a true_ide vector) so 
# every call site reads its own intent directly

power_rate <- function(p_vec) return(100 * mean(p_vec <= 0.05, na.rm = TRUE))

# mlm_exclusions(): because reject_rate() & power_rate() drop NA p-values from
# the denominator, this counts, for every p_* column of an MLM results frame,
# the replications with no usable p-value (is.na(p), i.e. excluded from the
# rate) & those whose conv_* flag is FALSE (lavaan reported non-convergence of
# either nested model; such a fit can still return a p-value, so the two counts
# need not agree). One row per evaluation; source names the results frame

mlm_exclusions <- function(res, source) {
  p_cols     <- grep("^p_", names(res), value = TRUE)
  key        <- sub("^p_", "", p_cols)
  conv_cols  <- paste0("conv_", key)
  out        <- data.frame(
    Source                = source,
    Scenario              = substr(key, 1, 1),
    Evaluation            = ifelse(grepl("vs_comparator", key), 
                                   "With Comparator", "Intervention Arm"),
    IDE                   = ifelse(grepl("_neg$", key), "-ve", "None"),
    Replications          = nrow(res),
    `No usable p-value`   = vapply(p_cols,    
                                   function(cl) return(sum(is.na(res[[cl]]))), 
                                   integer(1)),
    `Not converged`       = vapply(conv_cols, 
                                   function(cl) return(sum(!res[[cl]])),
                                   integer(1)),
    check.names = FALSE, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  return(out) }

## FUNCTIONS  -- LAVAAN MULTI-GROUP GROWTH MODEL (TWO-ARM, "MLM" LABEL) --------

# growth() fits latent intercept (i) & slope (s) factors for the 3 occasions,
# with factor loadings fixed to elapsed time in years - MUST BE CENTRED as per
# theory: (-1.5, 0, 1.5), so the slope factor is in yearly units & the
# intercept factor is anchored at the same middle occasion (age 2.5) throughout 
# the script. Cov(intercept, slope) is NOT invariant to where the intercept is 
# anchored, so an uncentred set of loadings (e.g. 0, 1.5, 3, anchoring at 
# baseline) tests a DIFFERENT covariance - centring at age 2.5 is necessary. 

# group = "arm" fits the model separately per arm; lavTestLRT() contrasts the 
# freely-estimated model (mod_free, intercept-slope covariance unconstrained per
# arm) against mod_equal (covariance constrained equal across arms via 
# group.equal). converged is TRUE only if BOTH models converged.

# RESIDUAL VARIANCE TIED WITHIN ARM (labels e1/e2): a free residual variance
# per occasion (growth()'s default) makes the per-arm covariance structure
# saturated, which under this population's lag-decay correlation structure
# forces a Heywood case (improper, usually negative wt_age1 residual
# variance) & invalidates the chi-square LRT. Per-group labels c(e1, e2) tie
# wt_age1/wt_age25/wt_age4's residual variances within each arm (one shared
# occasion-level residual per arm) but leave them free BETWEEN arms.
# se = "none" & baseline = FALSE (here & in lavaan_one_arm()) skip standard
# errors & the baseline (independence) model, neither of which the LRT uses:
# the p-value is unchanged & each fit runs about a third faster.

lavaan_two_arm <- function(pop) {
  pop$arm      <- factor(pop$group, levels = c(0, 1), 
                         labels = c("control", "treat") )
  model_syntax <- "
    i =~ 1*wt_age1 + 1*wt_age25 + 1*wt_age4
    s =~ -1.5*wt_age1 + 0*wt_age25 + 1.5*wt_age4
    wt_age1  ~~ c(e1, e2)*wt_age1
    wt_age25 ~~ c(e1, e2)*wt_age25
    wt_age4  ~~ c(e1, e2)*wt_age4 "
  fit_result   <- tryCatch({
    mod_free   <- growth(model_syntax, data = pop, group = "arm", se = "none", 
                         baseline = FALSE)
    mod_equal  <- growth(model_syntax, data = pop, group = "arm", 
                         group.equal = "lv.covariances", se = "none", 
                         baseline = FALSE)
    lr_test    <- lavTestLRT(mod_equal, mod_free)
    conv_free  <- lavInspect(mod_free,  "converged")
    conv_equal <- lavInspect(mod_equal, "converged")
    list(p_val = lr_test[2, "Pr(>Chisq)"], ok = TRUE, conv_free = conv_free, 
         conv_equal = conv_equal) },
    error = function(e) list(p_val = NA_real_, ok = FALSE, conv_free = FALSE, 
                             conv_equal = FALSE, err = conditionMessage(e) ) )
  converged    <- fit_result$ok && isTRUE(fit_result$conv_free) && 
    isTRUE(fit_result$conv_equal)
  return(list(p = fit_result$p_val, converged = converged, 
              conv_uncon = fit_result$conv_free, 
              conv_con = fit_result$conv_equal,
              err = fit_result$err) ) }

## FUNCTIONS  -- LAVAAN SINGLE-GROUP GROWTH MODEL (ONE-ARM, "MLM" LABEL) -------

# Lavaan analogue of run_oldham()'s intervention_only test. Oldham's test
# correlates (baseline - follow-up) with their average WITHIN ONE ARM, testing
# whether that correlation differs from zero - a direct test of unequal
# baseline/follow-up variance for that arm alone. This fits the SAME
# single-group growth model as lavaan_two_arm() (identical loadings, identical
# tied-residual structure) on one arm's data only, & tests whether the
# intercept-slope covariance differs from zero via LRT: free model vs a model
# with i ~~ 0*s (that one covariance fixed to zero) - the free-vs-constrained
# LRT pattern mirrors lavaan_two_arm()'s free-vs-equal contrast, just on one
# arm & one parameter instead of a cross-arm equality.

lavaan_one_arm <- function(pop, group_value = 1) {
  arm_pop      <- pop[pop$group == group_value, ]
  model_free   <- "
    i =~ 1*wt_age1 + 1*wt_age25 + 1*wt_age4
    s =~ -1.5*wt_age1 + 0*wt_age25 + 1.5*wt_age4
    wt_age1  ~~ e*wt_age1
    wt_age25 ~~ e*wt_age25
    wt_age4  ~~ e*wt_age4 "
  model_con    <- paste0(model_free, "\n    i ~~ 0*s ")
  fit_result   <- tryCatch({
    mod_free   <- growth(model_free, data = arm_pop, se = "none", 
                         baseline = FALSE)
    mod_con    <- growth(model_con,  data = arm_pop, se = "none", 
                         baseline = FALSE)
    lr_test    <- lavTestLRT(mod_con, mod_free)
    conv_free  <- lavInspect(mod_free, "converged")
    conv_con   <- lavInspect(mod_con,  "converged")
    list(p_val = lr_test[2, "Pr(>Chisq)"], ok = TRUE, conv_free = conv_free, 
         conv_con = conv_con) },
    error = function(e) list(p_val = NA_real_, ok = FALSE, conv_free = FALSE, 
                             conv_con = FALSE, err = conditionMessage(e) ) )
  converged    <- fit_result$ok && isTRUE(fit_result$conv_free) && 
    isTRUE(fit_result$conv_con)
  return(list(p = fit_result$p_val, converged = converged, 
              conv_uncon = fit_result$conv_free, 
              conv_con = fit_result$conv_con,
              err = fit_result$err) ) }

## FUNCTIONS  -- COMBINED CONTRASTS (SAME SHAPE AS run_oldham) -----------------

# Runs the "MLM" contrast (lavaan-powered throughout) on one allocated & 
# IDE-injected population, returning the same intervention_only / vs_comparator 
# shape as run_oldham(), so both methods slot into the driver Table 2 identically.

mlm_contrasts <- function(pop) {
  intervention_only <- lavaan_one_arm(pop, group_value = 1)
  contrast          <- lavaan_two_arm(pop)
  return(list(intervention_only = intervention_only, vs_comparator = contrast)) }

## FUNCTIONS  -- PER-REPLICATE FISHER Z CAPTURE (CELLS C & D) ------------------

# Same construction as run_one_replication_CD, but returns the raw Fisher z 
# values for both arms rather than collapsing to a p-value; this is the only 
# additional information needed for the diagnostic, & the per-replicate seeds 
# are identical to the sweep so the two are directly comparable 

run_one_replication_CD_z <- function(rep_i, z_tag) {
  pool_set       <- diag_pools[[z_tag]]
  draw_treat_C   <- sample_from_pool(pool_set$treat_C,   
                                     n_draw = Nsmp, 
                                     seed = seed_for("sweep", rep_i, "sample",
                                                     "treat",   "C"))
  draw_control_C <- sample_from_pool(pool_set$control_C, 
                                     n_draw = Nsmp, 
                                     seed = seed_for("sweep", rep_i,
                                                     "sample", "control", "C"))
  draw_treat_D   <- sample_from_pool(pool_set$treat_D,   
                                     n_draw = Nsmp, 
                                     seed = seed_for("sweep", rep_i, "sample", 
                                                     "treat",   "D"))
  draw_control_D <- sample_from_pool(pool_set$control_D, 
                                     n_draw = Nsmp, 
                                     seed = seed_for("sweep", rep_i, "sample", 
                                                     "control", "D"))
  old_treat_C    <- oldham_one_group(draw_treat_C$wt_age1,   
                                     draw_treat_C$wt_age4)
  old_control_C  <- oldham_one_group(draw_control_C$wt_age1, 
                                     draw_control_C$wt_age4)
  old_treat_D    <- oldham_one_group(draw_treat_D$wt_age1,   
                                     draw_treat_D$wt_age4)
  old_control_D  <- oldham_one_group(draw_control_D$wt_age1, 
                                     draw_control_D$wt_age4)
  return(data.frame(
    rep  = rep_i,
    z1_C = old_treat_C$z,
    z0_C = old_control_C$z,
    n1_C = old_treat_C$n,
    n0_C = old_control_C$n,
    z1_D = old_treat_D$z,
    z0_D = old_control_D$z,
    n1_D = old_treat_D$n,
    n0_D = old_control_D$n)) }

## FUNCTIONS  -- EMPIRICAL VS. NOMINAL VARIANCE OF FISHER Z CONTRAST -----------

# For one cell (C or D) at one threshold, computes empirical variance of (z1-z0) 
# across replicates & the nominal variance implied by 1/(n1-3) + 1/(n0-3) where 
# n1 & n0 are constant across replicates (= Nsmp) so the nominal variance is a 
# single number, not a vector

diagnose_cell <- function(z1, z0, n1, n0) {
  contrast             <- z1 - z0
  nominal_variance     <- 1 / (n1[1] - 3) + 1 / (n0[1] - 3)
  empirical_variance   <- var(contrast)
  return(list(
    empirical_variance = empirical_variance,
    nominal_variance   = nominal_variance,
    ratio              = empirical_variance / nominal_variance,
    mean_contrast      = mean(contrast))) }

## FUNCTIONS  -- PLOT SAVE & TABLE FORMAT HELPERS ------------------------------

# Plot save function (plot_w = 27 cm, plot_h = 17 cm) in Plots

SaveJpg <- function(plt, name, width = plot_w, height = plot_h, dpi = 600) {
  ggsave(plt, file = file.path(Plots, paste0(name, ".jpg")), device = "jpeg",
         width = width, height = height, dpi = dpi, limitsize = FALSE)
  return(invisible(NULL)) }

# Table 2 helpers. Table 2 (PART SIX) is long format: one row per cell x
# IDE status (None = no_ide, -ve = true_ide), rather than one row per cell x
# test with Type 1 error & power as side-by-side columns. Every row carries
# all four value columns (Oldham T1 Error, Oldham Power, MLM T1 Error, MLM
# Power), but only the column matching that row's IDE status is populated -
# the other three show "-" (not NA/blank), since a Type 1 error is only
# defined on a no_ide draw & a power figure is only defined on a true_ide
# draw; showing "-" makes that explicit rather than leaving an ambiguous
# empty cell. The Intervention Arm block (cells A-D) populates BOTH Oldham's
# single-arm test & its MLM analogue (lavaan_one_arm(), PART TWO - free vs
# covariance-constrained single-group growth model, the MLM equivalent of
# Oldham's within-arm baseline/follow-up correlation test). The With
# Comparator block (cells B-D only; no comparator arm exists for Cell A)
# populates both Oldham & MLM's two-arm test.

# pct_or_dash(): formats an already-computed rate as "XX.X%", or returns "-"
# when show is FALSE, so the same helper handles every value cell in the
# table without a long literal vector per column

pct_or_dash <- function(rate, show) if (show) sprintf("%.1f%%", rate) else "-"

# variance_label()/recruitment_label(): pulled from cell_grid (PART ONE) so
# the table's wording stays a single source of truth with the rest of the
# script, rather than a second hard-coded copy of the same regime/allocation
# mapping

variance_label    <- function(cell_id) 
  if (cell_grid$variance_regime[cell_grid$cell == cell_id] == "homo") "Constant" else "Non-constant"
recruitment_label <- function(cell_id) 
  if (cell_grid$allocation_rule[cell_grid$cell == cell_id] == "random") "Random" else "Targeted"

# table2_row(): one row of Table2 for a given cell, IDE status, & evaluation
# block (called from PART SIX, after verify_res & mlm_res exist).
# oldham_t1/oldham_pow/mlm_t1/mlm_pow are the already-computed rates (left NA
# where not applicable for that row); is_none controls which pair of columns
# (T1 Error vs Power) gets the real value & which gets "-"

table2_row <- function(evaluation, cell_id, ide_label, is_none, oldham_t1 = NA, 
                       oldham_pow = NA, mlm_t1 = NA, mlm_pow = NA) {
  data.frame(
    Evaluation        = evaluation,
    Cell              = cell_id,
    Variance          = variance_label(cell_id),
    Recruitment       = recruitment_label(cell_id),
    IDE               = ide_label,
    `Oldham T1 Error` = pct_or_dash(oldham_t1,  is_none && !is.na(oldham_t1) ),
    `Oldham Power`    = pct_or_dash(oldham_pow, !is_none && !is.na(oldham_pow) ),
    `MLM T1 Error`    = pct_or_dash(mlm_t1,     is_none && !is.na(mlm_t1) ),
    `MLM Power`       = pct_or_dash(mlm_pow,    !is_none && !is.na(mlm_pow) ),
    check.names       = FALSE, stringsAsFactors = FALSE) }

# tableA1_row(): one row of Table A1 (Appendix, alternative tracking-
# preserving true_ide) - power-only, since Type 1 error is identical to
# Table 2 under both true_ide specifications & is not repeated here.
# is_none = TRUE (the no_ide row) always shows "-" for both power columns:
# no_ide has no power to report, & that row exists only so the table shows
# the same IDE = no_ide/-ve row pairing as Table 2

tableA1_row <- function(evaluation, cell_id, ide_label, is_none,
                        oldham_pow = NA, mlm_pow = NA) {
  data.frame(
    Evaluation     = evaluation,
    Cell           = cell_id,
    Variance       = variance_label(cell_id),
    Recruitment    = recruitment_label(cell_id),
    IDE            = ide_label,
    `Oldham Power` = pct_or_dash(oldham_pow, !is_none && !is.na(oldham_pow) ),
    `MLM Power`    = pct_or_dash(mlm_pow,    !is_none && !is.na(mlm_pow) ),
    check.names    = FALSE, stringsAsFactors = FALSE) }

## #############################################################################
## PART THREE -- FUNCTIONS VERIFICATIONS
## PROGRAM    -- CONFIRM CELL GRID ---------------------------------------------

# This part is distinct from PART TWO's function definitions above & PART FOUR's 
# simulation/cache code below, so every executed PROGRAM step, as opposed to a 
# function definition, lives under its own PART, matching the FUNCTIONS vs. 
# PROGRAM split used throughout.

log_step("IDE PAPER: 2x2 CELL GRID:", 1)
print(cell_grid[, c("cell", "variance_regime", "allocation_rule", "simulate")])
log_step("Setup & global parameters loaded")

## PROGRAM    -- VERIFY POPULATION GENERATOR -----------------------------------

log_step("POPULATION GENERATOR: VERIFICATION", 1)

# Draw one population under each variance regime & print verification against
# CDC/LMS targets; confirms this section's own population generator (PART TWO's
# generate_population() & verify_population()) is correct standalone, BEFORE
# the more complex allocation/IDE logic tested later (PART FOUR's "VERIFY
# ALLOCATION & IDE SWITCHES", which runs after the targeted pools it needs are
# built) is layered on top - if this check fails, the fault is in population
# generation itself, not in allocation or IDE effects

mu_pop         <- c(mu_age1, mu_age25, mu_age4)
pop_homo_chk   <- generate_population(N = Nsmp, mu_vec = mu_pop, 
                                      sd_vec = sd_homo, rho_base = rho_base, 
                                      seed = seed_for("verify", "homo"))
pop_hetero_chk <- generate_population(N = Nsmp, mu_vec = mu_pop, 
                                      sd_vec = sd_hetero, 
                                      rho_base = rho_base, 
                                      seed = seed_for("verify", "hetero"))
log_step("homo regime:", 1)
print(verify_population(pop_homo_chk, mu_pop, sd_homo, rho_base)$moments)
print(verify_population(pop_homo_chk, mu_pop, sd_homo, rho_base)$correlations)
log_step("hetero regime:", 1)
print(verify_population(pop_hetero_chk, mu_pop, sd_hetero, rho_base)$moments)
print(verify_population(pop_hetero_chk, mu_pop, sd_hetero, rho_base)$correlations)

## #############################################################################
## PART FOUR  -- SIMULATION, CACHE, & DIAGNOSTICS
## PROGRAM    -- BUILD TARGETED-ALLOCATION POOLS (CELLS C & D) -----------------

log_step("TARGETED SELECTION POOLS: BUILD (CELLS C & D)", 1)

# Each targeted "pool" is the specification of one arm's baseline-truncated
# distribution (build_targeted_pool(), PART TWO), from which every
# replication draws its children afresh (sample_from_pool()); building one
# costs nothing, so nothing here is cached. Control pools are shared across
# ide_types (control arm is never touched by any intervention effect);
# treatment pools built separately PER ide_type, since the intervention effect 
# is baked directly into the generating distribution (mean shift for no_ide; 
# mean shift & weakened tracking for true_ide), not applied afterwards.

# One row per pool to build: cell, role (treat/control), variance regime,
# ide_type ("" for control, since control has none)

pool_spec <- data.frame(
  cell      = c("C",    "C",    "C",       "D",    "D",    "D"),
  role      = c("control", "treat", "treat", "control", "treat", "treat"),
  sd_regime = c("homo",    "homo",   "homo",   "hetero", "hetero", "hetero"),
  ide_type  = c("",     "no_ide", "true_ide",  "",     "no_ide", "true_ide"),
  stringsAsFactors = FALSE)

pools        <- list()
pool_tag     <- cache_tag(cache_version, mu_pop, sd_homo, sd_hetero, rho_base, 
                          mean_mult, sd_mult_true_ide)
for (row_i in seq_len(nrow(pool_spec))) {
  spec       <- pool_spec[row_i, ]
  sd_vec     <- if (spec$sd_regime == "hetero") sd_hetero else sd_homo
  pool_name  <- if (spec$role == "control") 
    sprintf("pool_control_%s", spec$cell) else
    sprintf("pool_treat_%s_%s", spec$cell, spec$ide_type)
  if (spec$role == "control") {
    pool     <- build_targeted_pool(mu_vec = mu_pop, sd_vec = sd_vec, 
                                    rho_base = rho_base, ref_mean = mu_pop, 
                                    ref_sd = sd_vec)
  } else {
    dist     <- build_treatment_distribution(mu_vec = mu_pop, sd_vec = sd_vec, 
                                             rho_base = rho_base, 
                                             mean_mult = mean_mult, 
                                             sd_mult_true_ide = sd_mult_true_ide, 
                                             ide_type = spec$ide_type)
    pool     <- build_targeted_pool(mu_vec = dist$mu_vec, sd_vec = sd_vec, 
                                    rho_base = rho_base, ref_mean = mu_pop, 
                                    ref_sd = sd_vec,
                                    Sigma_override = dist$Sigma) }
  pools[[pool_name]] <- pool }

pool_control_C       <- pools[["pool_control_C"]]
pool_control_D       <- pools[["pool_control_D"]]
pool_treat_C_no_ide  <- pools[["pool_treat_C_no_ide"]]
pool_treat_C_true_ide <- pools[["pool_treat_C_true_ide"]]
pool_treat_D_no_ide  <- pools[["pool_treat_D_no_ide"]]
pool_treat_D_true_ide <- pools[["pool_treat_D_true_ide"]]

log_step(sprintf("Baseline threshold (wt_age1 above): C = %.3f kg   D = %.3f kg",
                 pool_control_C$threshold, pool_control_D$threshold))

## PROGRAM    -- APPENDIX: ALTERNATIVE true_ide TREATMENT POOLS (CELLS C & D) --

log_step("APPENDIX: ALTERNATIVE (TRACKING-PRESERVING) true_ide POOLS (CELLS C & D)", 1)

# Only the two true_ide TREATMENT pools differ under the alternative
# (build_treat_dist_condshrink()) definition - control pools are
# identical either way (control is never touched by ide_type) & no_ide
# treatment pools are also identical either way (both builders return the
# SAME Sigma when ide_type == "no_ide"), so pool_control_C/D & pool_treat_C/D
# _no_ide above are reused unchanged; only pool_treat_C_true_ide & pool_treat_
# D_true_ide need alternative versions, built here for the Appendix Table 2
# robustness check (PART TWO's build_treat_dist_condshrink()).

pool_tag_alt <- cache_tag(cache_version, mu_pop, sd_homo, sd_hetero, rho_base,
                          mean_mult, sd_mult_true_ide, "condshrink")
if (!skip_appendix) {
  appendix_pool_spec <- pool_spec[pool_spec$ide_type == "true_ide", ]
  appendix_pools     <- list()
  for (row_i in seq_len(nrow(appendix_pool_spec))) {
    spec       <- appendix_pool_spec[row_i, ]
    sd_vec     <- if (spec$sd_regime == "hetero") sd_hetero else sd_homo
    pool_name  <- sprintf("pool_treat_%s_%s", spec$cell, spec$ide_type)
    dist       <- build_treat_dist_condshrink(mu_vec = mu_pop, sd_vec = sd_vec,
                                              rho_base = rho_base, 
                                              mean_mult = mean_mult,
                                              sd_mult_true_ide = sd_mult_true_ide,
                                              ide_type = spec$ide_type)
    appendix_pools[[pool_name]] <- build_targeted_pool(mu_vec = dist$mu_vec, 
                                                       sd_vec = sd_vec,
                                                       rho_base = rho_base, 
                                                       ref_mean = mu_pop, 
                                                       ref_sd = sd_vec,
                                                       Sigma_override = dist$Sigma) }
  pool_treat_C_true_ide_alt <- appendix_pools[["pool_treat_C_true_ide"]]
  pool_treat_D_true_ide_alt <- appendix_pools[["pool_treat_D_true_ide"]]
} else { log_step("skip_appendix = TRUE -- skipping alternative true_ide pools (C & D)") }

## PROGRAM    -- VERIFY ALLOCATION & IDE SWITCHES ------------------------------

log_step("ALLOCATION & IDE SWITCHES: VERIFICATION", 1)

# Runs after the targeted pools above are built, since Cell C's demonstration
# draw below reads pool_treat_C_no_ide / pool_control_C directly.

# Cell B (random/universal): TWO independent N = 500 draws from the SAME
# heteroscedastic population parametrisation -- one becomes the treatment
# arm, one becomes the control arm. Not a single population split by coin
# flip; the two arms are different children, drawn separately, exactly as
# for the targeted cells (C, D), but without any pool/threshold step since
# allocation here is universal/untargeted. The treatment arm's distribution
# (mean shift, & for true_ide, weakened tracking) is baked directly into
# the generating Sigma via build_treatment_distribution(), not applied
# afterward as a reconstruction

dist_B_no_ide  <- build_treatment_distribution(mu_vec = mu_pop, 
                                               sd_vec = sd_hetero, 
                                               rho_base = rho_base,
                                               mean_mult = mean_mult, 
                                               sd_mult_true_ide = sd_mult_true_ide,
                                               ide_type = "no_ide")
draw_treat_B   <- generate_population(N = Nsmp, mu_vec = dist_B_no_ide$mu_vec, 
                                      sd_vec = sd_hetero, rho_base = rho_base, 
                                      seed = seed_for("pop", "hetero", "treat"),
                                      Sigma_override = dist_B_no_ide$Sigma)
draw_control_B <- generate_population(N = Nsmp, mu_vec = mu_pop, 
                                      sd_vec = sd_hetero, rho_base = rho_base, 
                                      seed = seed_for("pop", "hetero", "control"))
draw_treat_B$group   <- 1
draw_control_B$group <- 0
cell_B <- rbind(draw_control_B, draw_treat_B)
log_step(sprintf("Cell B - random   draw (control = %d, treat = %d) - two independent N = 500 draws",
                 sum(cell_B$group == 0), sum(cell_B$group == 1)))

# Cell C (targeted): group0 & group1 are TWO independent N = 500 draws from its
# own at-risk truncated pool (pool_treat_C_no_ide, pool_control_C); both pools
# are thresholded IDENTICALLY (top ~5% of baseline weight, z > 1.645), since an
# exchangeable comparator must share the treatment arm's restricted baseline
# variance structure. The treatment pool already has the no_ide intervention
# effect baked into its generating distribution. The two arms are drawn
# independently, which reflects that the treatment
# & control arms are different children

draw_treat_C   <- sample_from_pool(pool_treat_C_no_ide, 
                                   n_draw = Nsmp, 
                                   seed = seed_for("sample", "treat",   "C"))
draw_control_C <- sample_from_pool(pool_control_C,      
                                   n_draw = Nsmp, 
                                   seed = seed_for("sample", "control", "C"))
draw_treat_C$group   <- 1
draw_control_C$group <- 0
cell_C <- rbind(draw_control_C, draw_treat_C)
log_step(sprintf("Cell C - targeted draw (control = %d, treat = %d) - two independent N = 500 pool draws",
                 sum(cell_C$group == 0), sum(cell_C$group == 1)))
log_step(sprintf("Cell C - mean wt_age1 group0 = %.3f  group1 = %.3f  (both above at-risk threshold, drawn independently)",
                 mean(cell_C$wt_age1[cell_C$group == 0]),
                 mean(cell_C$wt_age1[cell_C$group == 1])))

# Confirm Cell B's true_ide treatment pool has a smaller follow-up SD than its
# no_ide counterpart while the control arm (untouched by any ide_type) is
# identical regardless. Both treatment distributions are regenerated fresh
# here for the demonstration, since the effect is baked in at generation
# time rather than modifiable after the fact

dist_B_true_ide  <- build_treatment_distribution(mu_vec = mu_pop, 
                                                 sd_vec = sd_hetero, 
                                                 rho_base = rho_base,
                                                 mean_mult = mean_mult, 
                                                 sd_mult_true_ide = sd_mult_true_ide,
                                                 ide_type = "true_ide")
draw_treat_B_neg <- generate_population(N = Nsmp, 
                                        mu_vec = dist_B_true_ide$mu_vec, 
                                        sd_vec = sd_hetero,
                                        rho_base = rho_base, 
                                        seed = seed_for("pop", "hetero", 
                                                        "treat", "true_ide"),
                                        Sigma_override = dist_B_true_ide$Sigma)
log_step(sprintf("no_ide  intervention-arm SD at age4 = %.3f", 
                 sd(draw_treat_B$wt_age4)))
log_step(sprintf("true_ide intervention-arm SD at age4 = %.3f   (should be smaller than no_ide)", 
                 sd(draw_treat_B_neg$wt_age4)))
log_step(sprintf("control-arm SD at age4 identical across ide_type:   no_ide = %.3f   true_ide = %.3f", 
                 sd(draw_control_B$wt_age4), sd(draw_control_B$wt_age4)))

## PROGRAM    -- OLDHAM MODULE (K REPLICATIONS, ALL FOUR CELLS) ----------------

log_step(sprintf("OLDHAM MODULE (K = %d REPLICATIONS, ALL FOUR CELLS) ", K), 1)

# K is the number of replications. A single draw cannot distinguish a genuine
# finding from sampling noise for anything involving a p-value or rejection 
# rate, so calibration checks loop K times & report rejection proportions.
# This dispatch must run AFTER the targeted pools above are built, since
# run_one_replication()'s default pool_set argument (PART TWO) reads
# pool_treat_C_no_ide etc. directly from this environment.

# Treatment-arm with no IDE for A & B (no_ide): fixed across all K replicates,
# computed once here rather than inside the per-replicate function, since
# recomputing an identical Sigma/mu_int K times would be wasted work

dist_A_no_ide <- build_treatment_distribution(mu_vec = mu_pop, sd_vec = sd_homo, 
                                              rho_base = rho_base,
                                              mean_mult = mean_mult, 
                                              sd_mult_true_ide = sd_mult_true_ide,
                                              ide_type = "no_ide")
dist_B_no_ide <- build_treatment_distribution(mu_vec = mu_pop, 
                                              sd_vec = sd_hetero, 
                                              rho_base = rho_base,
                                              mean_mult = mean_mult, 
                                              sd_mult_true_ide = sd_mult_true_ide,
                                              ide_type = "no_ide")

# Treatment-arm with IDE for A & B (true_ide): the power counterpart of
# dist_A_no_ide/dist_B_no_ide above, across all K replicates

dist_A_true_ide <- build_treatment_distribution(mu_vec = mu_pop, 
                                                sd_vec = sd_homo, 
                                                rho_base = rho_base,
                                                mean_mult = mean_mult, 
                                                sd_mult_true_ide = sd_mult_true_ide,
                                                ide_type = "true_ide")
dist_B_true_ide <- build_treatment_distribution(mu_vec = mu_pop, 
                                                sd_vec = sd_hetero, 
                                                rho_base = rho_base,
                                                mean_mult = mean_mult, 
                                                sd_mult_true_ide = sd_mult_true_ide,
                                                ide_type = "true_ide")

# APPENDIX ONLY: cells A & B's true_ide distribution under the alternative
# (tracking-preserving) definition - see build_treatment_distribution_
# condshrink() (PART TWO) for what this represents. no_ide is unaffected by
# which builder is used, so dist_A_no_ide/dist_B_no_ide above are reused
# unchanged for the appendix run - only the true_ide distributions differ.

if (!skip_appendix) {
  dist_A_true_ide_alt <- build_treat_dist_condshrink(mu_vec = mu_pop, 
                                                     sd_vec = sd_homo, 
                                                     rho_base = rho_base,
                                                     mean_mult = mean_mult, 
                                                     sd_mult_true_ide = sd_mult_true_ide,
                                                     ide_type = "true_ide")
  dist_B_true_ide_alt <- build_treat_dist_condshrink(mu_vec = mu_pop, 
                                                     sd_vec = sd_hetero, 
                                                     rho_base = rho_base,
                                                     mean_mult = mean_mult, 
                                                     sd_mult_true_ide = sd_mult_true_ide,
                                                     ide_type = "true_ide") }
verify_tag       <- cache_tag(K, Nsmp, pool_tag, dist_A_no_ide$mu_vec, 
                              dist_B_no_ide$mu_vec, sd_homo, sd_hetero, 
                              rho_base, mu_pop)
verify_res_file  <- file.path(Data, sprintf("verify results %s.rds", verify_tag))
t_verify         <- Sys.time()
if (force_rerun_OLDHAM || !file.exists(verify_res_file)) {
  log_step(sprintf("no cache found (or force_rerun_OLDHAM = TRUE) -- running %d replications now", K))
  verify_res     <- run_reps_chunked(K, parallel_K, 
                                             chunk_size = chunk_size, 
                                             cl = get_script_cl(),
                                             checkpoint_tag = sprintf("verify %s",
                                                                      verify_tag),
                                             force_fresh = force_rerun_OLDHAM)
  saveRDS(verify_res, verify_res_file) 
  log_done(t_verify, sprintf("K = %d replications", K)) 
} else { 
  log_step(sprintf("loading cached results") )
  verify_res  <- readRDS(verify_res_file) }

# Cell-by-cell report of the K-replication run above (Type 1 error, power),
# using verify_res -- confirmed correctly positioned after PART TWO's function
# definitions, PART THREE's population/allocation/IDE checks, & this section's
# own pool-dependent dispatch, all of which verify_res depends on.

# Cell A: positive control. Intervention-only rejection rate should sit near 5%
# (correlation centred on zero, since baseline & follow-up SDs are equal) even 
# with no comparator; power against a genuine negative IDE is read from the same 
# intervention-only test applied to the true_ide draw

# Cell B: Intervention-only rejection rate should be close to 100% since 
# baseline & follow-up variances differ irrespective of any IDE; the comparator 
# contrast should restore Type 1 error close to 5%

# Cell C: homoscedastic + targeted, the extreme replication of Beggs et al. Even 
# with constant variance, truncated baseline allocation compresses each arm's 
# baseline variance more than its follow-up variance follow-up correlates 
# imperfectly with baseline), so intervention-only rejects near 100% just as 
# under heteroscedasticity. 

# Cell D: heteroscedastic + targeted, compounded case combining both challenges. 
# Intervention-only rejects near 100% for the same reasons as B & C.

cat(sprintf("           Cell A:
           homoscedastic, random, intervention-only     -- Type 1 error = %.1f%%    (expected ~5%%)",
            reject_rate(verify_res$p_A_int_only)), "\n")
cat(sprintf("           homoscedastic, random, intervention-only     -- Power        = %.1f%%",
            power_rate(verify_res$p_A_int_only_neg)), "\n")
cat(sprintf("           Cell B:
           heteroscedastic, random, intervention-only   -- Type 1 error = %.1f%%  (expected ~100%%)",
            reject_rate(verify_res$p_B_int_only)), "\n")
cat(sprintf("           heteroscedastic, random, intervention-only   -- Power        = %.1f%%",
            power_rate(verify_res$p_B_int_only_neg)), "\n")
cat(sprintf("           Cell C:
           homoscedastic, targeted, intervention-only   -- Type 1 error = %.1f%%  (expected ~100%%)",
            reject_rate(verify_res$p_C_int_only)), "\n")
cat(sprintf("           homoscedastic, targeted, intervention-only   -- Power        = %.1f%%",
            power_rate(verify_res$p_C_int_only_neg)), "\n")
cat(sprintf("           Cell D:
           heteroscedastic, targeted, intervention-only -- Type 1 error = %.1f%%  (expected ~100%%)",
            reject_rate(verify_res$p_D_int_only)), "\n")
cat(sprintf("           heteroscedastic, targeted, intervention-only -- Power        = %.1f%%",
            power_rate(verify_res$p_D_int_only_neg)), "\n\n")

cat(sprintf("           Cell B vs. comparator:
           heteroscedastic, random, with comparator     -- Type 1 error = %.1f%%    (expected ~5%%)",
            reject_rate(verify_res$p_B_vs_comparator)), "\n")
cat(sprintf("           heteroscedastic, random, with comparator     -- Power        = %.1f%%",
            power_rate(verify_res$p_B_vs_comparator_neg)), "\n")
cat(sprintf("           Cell C vs. comparator:
           homoscedastic, targeted, with comparator     -- Type 1 error = %.1f%%   (expected >5%% - biased)",
            reject_rate(verify_res$p_C_vs_comparator)), "\n")
cat(sprintf("           homoscedastic, targeted, with comparator     -- Power        = %.1f%%",
            power_rate(verify_res$p_C_vs_comparator_neg)), "\n")
cat(sprintf("           Cell D vs. comparator:
           heteroscedastic, targeted, with comparator   -- Type 1 error = %.1f%%   (expect biased, > 5%%)",
            reject_rate(verify_res$p_D_vs_comparator)), "\n")
cat(sprintf("           heteroscedastic, targeted, with comparator   -- Power        = %.1f%%",
            power_rate(verify_res$p_D_vs_comparator_neg)), "\n")

log_step("The Oldham module is verified across all four cells for Type 1 error & Power")

## PROGRAM    -- MLM MODULE (K REPLICATIONS, ALL FOUR CELLS) -------------------

log_step(sprintf("MLM MODULE (K = %d REPLICATIONS, ALL FOUR CELLS, PARALLEL)", K), 1)

# K is shared with the Oldham pass, the threshold sweep, & the diagnostic; MLM
# will still be slower than Oldham at any given K, since each lavaan growth-
# model fit is far more expensive per replication than Oldham's closed-form
# correlation, even run in parallel. One replication draws the identical
# children as the matching Oldham replication - same seed keys as
# run_one_replication() - so the 2 methods are evaluated on same data & are
# comparable cell by cell.

# MLM module verification: MLM replications across all four cells. Population
# draws are byte-for byte identical to run_one_replication() - same seed keys -
# so MLM & Oldham passes are evaluated on exactly the same children & are thus
# comparable. Cell A has no comparator arm at all (no group == 0 draw exists),
# so it goes through lavaan_one_arm() directly rather than mlm_contrasts()
# (which expects both arms present for its vs_comparator leg); cells B-D go
# through mlm_contrasts(), which populates BOTH intervention_only (via
# lavaan_one_arm() on group == 1) & vs_comparator (via lavaan_two_arm()).
# lavaan runs fully in-process (no external executable, no work_dir/log files
# to isolate between replications), so no mlm_work_root/failed-log machinery
# is needed here - each replication's tryCatch (inside lavaan_one_arm()/
# lavaan_two_arm()) already handles fit failures cleanly on its own.
# The random-allocation cells (A, B) & the targeted cells (C, D) run as two
# separately cached passes, merged by replication: A & B never touch the
# targeted draws, so a change to the targeted-recruitment machinery re-runs
# only the C & D pass. Each replication's seeds depend only on rep_i & the
# cell, so the split gives exactly the same results as one combined pass.
# dist_A_neg_arg / dist_B_neg_arg: see matching note on run_one_replication()
# above - defaults to the primary globals, overridable for the Appendix run.

run_one_replication_mlm_AB <- function(rep_i, dist_A_neg_arg = dist_A_true_ide, 
                                       dist_B_neg_arg = dist_B_true_ide) {
  pop_A     <- generate_population(N = Nsmp, mu_vec = dist_A_no_ide$mu_vec, 
                                   sd_vec = sd_homo, rho_base = rho_base, 
                                   seed = seed_for("verify", rep_i, "pop", "A"),
                                   Sigma_override = dist_A_no_ide$Sigma)
  pop_A_neg <- generate_population(N = Nsmp, mu_vec = dist_A_neg_arg$mu_vec, 
                                   sd_vec = sd_homo, 
                                   rho_base = rho_base, 
                                   seed = seed_for("verify", rep_i, "pop", "A", 
                                                   "neg"),
                                   Sigma_override = dist_A_neg_arg$Sigma)
  pop_A$group <- 1; pop_A_neg$group <- 1
  draw_treat_B     <- generate_population(N = Nsmp, 
                                          mu_vec = dist_B_no_ide$mu_vec, 
                                          sd_vec = sd_hetero, 
                                          rho_base = rho_base, 
                                          seed = seed_for("verify", rep_i,
                                                          "pop", "B", "treat"), 
                                          Sigma_override = dist_B_no_ide$Sigma)
  draw_treat_B_neg <- generate_population(N = Nsmp, 
                                          mu_vec = dist_B_neg_arg$mu_vec, 
                                          sd_vec = sd_hetero,
                                          rho_base = rho_base, 
                                          seed = seed_for("verify", rep_i, 
                                                          "pop", "B", "treat",
                                                          "neg"),
                                          Sigma_override = dist_B_neg_arg$Sigma)
  draw_control_B   <- generate_population(N = Nsmp, mu_vec = mu_pop, 
                                          sd_vec = sd_hetero,
                                          rho_base = rho_base, 
                                          seed = seed_for("verify", rep_i, 
                                                          "pop", "B", "control"))
  draw_treat_B$group <- 1; draw_control_B$group <- 0; draw_treat_B_neg$group <- 1
  cA      <- pop_A
  cA_neg  <- pop_A_neg
  cB      <- rbind(draw_control_B, draw_treat_B)
  cB_neg  <- rbind(draw_control_B, draw_treat_B_neg)
  mA      <- lavaan_one_arm(cA,     group_value = 1)
  mA_neg  <- lavaan_one_arm(cA_neg, group_value = 1)
  mB      <- mlm_contrasts(cB)
  mB_neg  <- mlm_contrasts(cB_neg)
  return(data.frame(
    rep                      = rep_i,
    p_A_int_only             = mA$p,
    conv_A_int_only          = mA$converged,
    p_A_int_only_neg         = mA_neg$p,
    conv_A_int_only_neg      = mA_neg$converged,
    p_B_int_only             = mB$intervention_only$p,
    conv_B_int_only          = mB$intervention_only$converged,
    p_B_int_only_neg         = mB_neg$intervention_only$p,
    conv_B_int_only_neg      = mB_neg$intervention_only$converged,
    p_B_vs_comparator        = mB$vs_comparator$p,
    conv_B_vs_comparator     = mB$vs_comparator$converged,
    p_B_vs_comparator_neg    = mB_neg$vs_comparator$p,
    conv_B_vs_comparator_neg = mB_neg$vs_comparator$converged) ) }

run_one_replication_mlm_CD <- function(rep_i, pool_set = list(
  treat_C     = pool_treat_C_no_ide,
  treat_C_neg = pool_treat_C_true_ide,
  control_C   = pool_control_C,
  treat_D     = pool_treat_D_no_ide,
  treat_D_neg = pool_treat_D_true_ide,
  control_D   = pool_control_D)) {
  draw_treat_C           <- sample_from_pool(pool_set$treat_C,     
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", rep_i, 
                                                             "sample", "treat",
                                                             "C") )
  draw_treat_C_neg       <- sample_from_pool(pool_set$treat_C_neg, 
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", rep_i,
                                                             "sample", "treat",
                                                             "C", "neg"))
  draw_control_C         <- sample_from_pool(pool_set$control_C,   
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", rep_i,
                                                             "sample",
                                                             "control", "C"))
  draw_treat_D           <- sample_from_pool(pool_set$treat_D,     
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", rep_i,
                                                             "sample", "treat",
                                                             "D"))
  draw_treat_D_neg       <- sample_from_pool(pool_set$treat_D_neg, 
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", rep_i,
                                                             "sample", "treat",
                                                             "D", "neg"))
  draw_control_D         <- sample_from_pool(pool_set$control_D,   
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", rep_i,
                                                             "sample",
                                                             "control", "D"))
  draw_treat_C$group     <- 1; draw_control_C$group <- 0
  draw_treat_C_neg$group <- 1
  draw_treat_D$group     <- 1; draw_control_D$group <- 0
  draw_treat_D_neg$group <- 1
  cC      <- rbind(draw_control_C, draw_treat_C)
  cC_neg  <- rbind(draw_control_C, draw_treat_C_neg)
  cD      <- rbind(draw_control_D, draw_treat_D)
  cD_neg  <- rbind(draw_control_D, draw_treat_D_neg)
  mC      <- mlm_contrasts(cC)
  mC_neg  <- mlm_contrasts(cC_neg)
  mD      <- mlm_contrasts(cD)
  mD_neg  <- mlm_contrasts(cD_neg)
  return(data.frame(
    rep                      = rep_i,
    p_C_int_only             = mC$intervention_only$p,
    conv_C_int_only          = mC$intervention_only$converged,
    p_C_int_only_neg         = mC_neg$intervention_only$p,
    conv_C_int_only_neg      = mC_neg$intervention_only$converged,
    p_C_vs_comparator        = mC$vs_comparator$p,
    conv_C_vs_comparator     = mC$vs_comparator$converged,
    p_C_vs_comparator_neg    = mC_neg$vs_comparator$p,
    conv_C_vs_comparator_neg = mC_neg$vs_comparator$converged,
    p_D_int_only             = mD$intervention_only$p,
    conv_D_int_only          = mD$intervention_only$converged,
    p_D_int_only_neg         = mD_neg$intervention_only$p,
    conv_D_int_only_neg      = mD_neg$intervention_only$converged,
    p_D_vs_comparator        = mD$vs_comparator$p,
    conv_D_vs_comparator     = mD$vs_comparator$converged,
    p_D_vs_comparator_neg    = mD_neg$vs_comparator$p,
    conv_D_vs_comparator_neg = mD_neg$vs_comparator$converged) ) }

# Power-only variants of run_one_replication_mlm_AB() / _CD() for the
# Appendix: same rationale as run_one_rep_power_only() above (PART
# TWO) - the Appendix's Type 1 error figures are reused from the primary
# mlm_res rather than recomputed, so only the true_ide lavaan fits (the
# expensive part) are run here, & only those Table A1 reports: Cell A's
# intervention-only test (lavaan_one_arm()) & the with-comparator contrast
# (lavaan_two_arm()) for cells B-D.

run_one_rep_mlm_power_only_AB <- function(rep_i, 
                                                  dist_A_neg_arg = dist_A_true_ide, 
                                                  dist_B_neg_arg = dist_B_true_ide) {
  pop_A_neg        <- generate_population(N = Nsmp, 
                                          mu_vec = dist_A_neg_arg$mu_vec, 
                                          sd_vec = sd_homo, rho_base = rho_base, 
                                          seed = seed_for("verify", rep_i, 
                                                          "pop", "A", "neg"),
                                          Sigma_override = dist_A_neg_arg$Sigma)
  pop_A_neg$group  <- 1
  draw_treat_B_neg <- generate_population(N = Nsmp, 
                                          mu_vec = dist_B_neg_arg$mu_vec, 
                                          sd_vec = sd_hetero, 
                                          rho_base = rho_base, 
                                          seed = seed_for("verify", rep_i,
                                                          "pop", "B", "treat",
                                                          "neg"),
                                          Sigma_override = dist_B_neg_arg$Sigma)
  draw_control_B   <- generate_population(N = Nsmp, mu_vec = mu_pop, 
                                          sd_vec = sd_hetero, 
                                          rho_base = rho_base, 
                                          seed = seed_for("verify", rep_i, "pop",
                                                          "B", "control"))
  draw_treat_B_neg$group <- 1; draw_control_B$group <- 0
  cA_neg  <- pop_A_neg
  cB_neg  <- rbind(draw_control_B, draw_treat_B_neg)
  mA_neg  <- lavaan_one_arm(cA_neg, group_value = 1)
  mB_neg  <- lavaan_two_arm(cB_neg)
  return(data.frame(
    rep                      = rep_i,
    p_A_int_only_neg         = mA_neg$p,
    conv_A_int_only_neg      = mA_neg$converged,
    p_B_vs_comparator_neg    = mB_neg$p,
    conv_B_vs_comparator_neg = mB_neg$converged) ) }

run_one_rep_mlm_power_only_CD <- function(rep_i, pool_set = list(
  treat_C_neg = pool_treat_C_true_ide,
  control_C   = pool_control_C,
  treat_D_neg = pool_treat_D_true_ide,
  control_D   = pool_control_D)) {
  draw_treat_C_neg       <- sample_from_pool(pool_set$treat_C_neg, 
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", rep_i, 
                                                             "sample", "treat",
                                                             "C", "neg"))
  draw_control_C         <- sample_from_pool(pool_set$control_C,   
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", rep_i,
                                                             "sample",
                                                             "control", "C"))
  draw_treat_D_neg       <- sample_from_pool(pool_set$treat_D_neg, 
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", rep_i,
                                                             "sample", "treat",
                                                             "D", "neg"))
  draw_control_D         <- sample_from_pool(pool_set$control_D,   
                                             n_draw = Nsmp, 
                                             seed = seed_for("verify", rep_i,
                                                             "sample",
                                                             "control", "D"))
  draw_treat_C_neg$group <- 1; draw_control_C$group <- 0
  draw_treat_D_neg$group <- 1; draw_control_D$group <- 0
  cC_neg  <- rbind(draw_control_C, draw_treat_C_neg)
  cD_neg  <- rbind(draw_control_D, draw_treat_D_neg)
  mC_neg  <- lavaan_two_arm(cC_neg)
  mD_neg  <- lavaan_two_arm(cD_neg)
  return(data.frame(
    rep                      = rep_i,
    p_C_vs_comparator_neg    = mC_neg$p,
    conv_C_vs_comparator_neg = mC_neg$converged,
    p_D_vs_comparator_neg    = mD_neg$p,
    conv_D_vs_comparator_neg = mD_neg$converged) ) }

# Runs one cached MLM pass: loads res_file if present (unless force_fresh),
# otherwise runs K replications of rep_fn in parallel, checkpointed under
# checkpoint_tag, & saves the result to res_file. MLM (lavaan) replications
# are slower per unit than Oldham/sweep/diagnostic, so global chunk_size
# (5,000) would leave a large amount of completed work exposed to loss if a
# crash happens mid-chunk. mlm_chnk_size gives MLM its own, smaller
# checkpoint granularity - a crash loses at most this many replications
# rather than up to chunk_size. MLM forces parallel_K = 1 so it always
# dispatches in parallel regardless of the top-level K; needs_parallel =
# TRUE tells get_script_cl() this call needs a real cluster even if K <=
# parallel_K at the top level (otherwise it would cache a NULL cluster
# decision on this, its first call, & every MLM-type stage would build/tear
# down its own private cluster instead of sharing one - see get_script_cl()'s
# note). Extra arguments (...) pass through to rep_fn.

run_mlm_pass <- function(res_file, rep_fn, checkpoint_tag, label, force_fresh, 
                         ...) {
  t_pass <- Sys.time()
  if (force_fresh || !file.exists(res_file) ) {
    log_step(sprintf("no cache found (or force rerun) -- running %d %s replications now (parallel)", 
                     K, label) )
    res  <- run_reps_chunked(K, parallel_K = 1,
                                     chunk_size = mlm_chnk_size,
                                     rep_fn = rep_fn, 
                                     cl = get_script_cl(needs_parallel = TRUE),
                                     checkpoint_tag = checkpoint_tag, 
                                     force_fresh = force_fresh, ...)
    saveRDS(res, res_file)
    log_done(t_pass, sprintf("K = %d %s replications", K, label) )
  } else {
    log_step(sprintf("loading cached %s results", label) )
    res     <- readRDS(res_file) }
  return(res) }

mlm_tag_AB <- cache_tag(cache_version, K, Nsmp, dist_B_no_ide$mu_vec, sd_homo, 
                        sd_hetero, rho_base, mu_pop,
                        "lavaan only v4", "cells A & B")
mlm_tag_CD <- cache_tag(cache_version, K, Nsmp, pool_tag, dist_B_no_ide$mu_vec,
                        sd_homo, sd_hetero, rho_base, mu_pop,
                        "lavaan only v4", "cells C & D")
mlm_res_AB <- run_mlm_pass(file.path(Data, 
                                     sprintf("mlm results %s.rds", mlm_tag_AB)), 
                           rep_fn = run_one_replication_mlm_AB,
                           checkpoint_tag = sprintf("mlm %s", mlm_tag_AB), 
                           label = "MLM (cells A & B)", 
                           force_fresh = force_rerun_MLM)
mlm_res_CD <- run_mlm_pass(file.path(Data, 
                                     sprintf("mlm results %s.rds", mlm_tag_CD)),
                           rep_fn = run_one_replication_mlm_CD,
                                    checkpoint_tag = sprintf("mlm %s", mlm_tag_CD),
                           label = "MLM (cells C & D)", 
                           force_fresh = force_rerun_MLM)
mlm_res    <- merge(mlm_res_AB, mlm_res_CD, by = "rep")
log_step("MLM module verified, cells A to D, Type 1 error & power (lavaan-powered growth-model LRT, PART TWO)")

## PROGRAM    -- APPENDIX: ALTERNATIVE true_ide ROBUSTNESS CHECK (TABLE A1 ONLY) --

if (skip_appendix) {
  log_step("skip_appendix = TRUE -- skipping Appendix module entirely (Table A1 will be omitted)", 1)
} else {
  log_step(sprintf("APPENDIX MODULE (K = %d REPLICATIONS, ALL FOUR CELLS): ALTERNATIVE (TRACKING-PRESERVING) true_ide",
                   K), 1)
  # Lighter-touch robustness check for the Appendix: reruns ONLY the POWER side
  # (true_ide) of the Table 2 comparison (Oldham + MLM, cells A-D) under the
  # alternative true_ide definition (build_treat_dist_condshrink(),
  # PART TWO) - not the threshold severity sweep or diagnostic, which are
  # specific to the targeted-recruitment story rather than the tracking-
  # mechanism assumption this appendix is stress-testing. Type 1 error is NOT
  # recomputed here at all: no_ide p-values are IDENTICAL to the primary run
  # regardless of which true_ide builder is used (both return the same Sigma for
  # ide_type == "no_ide"), so the Appendix Table 2 assembly below (PART SIX)
  # reuses verify_res/mlm_res's no_ide columns directly rather than paying for
  # a second, redundant no_ide pass - this was previously computed twice
  # (identically, up to Monte Carlo noise) & roughly doubled this section's
  # cost for no informational gain. Uses run_one_rep_power_only() &
  # run_one_rep_mlm_power_only_AB() / _CD(), which only draw & test the
  # true_ide arm, via their dist_A_neg_arg/dist_B_neg_arg override parameters,
  # with appendix_pool_set substituting the alternative-definition true_ide
  # pools for cells C & D built earlier in this PART (no_ide pools are no
  # longer needed here at all, so are dropped from this pool_set).
  appendix_pool_set <- list(
    treat_C_neg = pool_treat_C_true_ide_alt,
    control_C   = pool_control_C,
    treat_D_neg = pool_treat_D_true_ide_alt,
    control_D   = pool_control_D)
  verify_tag_alt      <- cache_tag(K, Nsmp, pool_tag_alt, dist_A_no_ide$mu_vec, 
                                   dist_B_no_ide$mu_vec,
                                   sd_homo, sd_hetero, rho_base, mu_pop,
                                   "condshrink", "power only")
  verify_res_alt_file <- file.path(Data, sprintf("verify results %s.rds", 
                                                 verify_tag_alt))
  t_verify_alt     <- Sys.time()
  if (force_rerun_OLDHAM || !file.exists(verify_res_alt_file)) {
    log_step(sprintf("no cache found (or force_rerun_OLDHAM = TRUE) -- running %d appendix Oldham power-only replications now", K))
    verify_res_alt <- run_reps_chunked(K, parallel_K, chunk_size = chunk_size, 
                                       cl = get_script_cl(),
                                       rep_fn         = run_one_rep_power_only,
                                       pool_set       = appendix_pool_set,
                                       dist_A_neg_arg = dist_A_true_ide_alt,
                                       dist_B_neg_arg = dist_B_true_ide_alt,
                                       checkpoint_tag = sprintf("verify alt %s",
                                                                verify_tag_alt),
                                       force_fresh    = force_rerun_OLDHAM)
    saveRDS(verify_res_alt, verify_res_alt_file)
    log_done(t_verify_alt, sprintf("K = %d appendix Oldham power-only replications", K))
  } else {
    log_step("loading cached appendix Oldham results")
    verify_res_alt <- readRDS(verify_res_alt_file) }
  # MLM power-only passes, split into cells A & B & cells C & D as for the
  # primary MLM pass above (run_mlm_pass())
  mlm_tag_alt_AB <- cache_tag(cache_version, K, Nsmp, dist_B_no_ide$mu_vec,
                              sd_homo, sd_hetero, rho_base, mu_pop,
                              "lavaan only v4 condshrink power only",
                              "cells A & B")
  mlm_tag_alt_CD <- cache_tag(cache_version, K, Nsmp, pool_tag_alt, 
                              dist_B_no_ide$mu_vec,
                              sd_homo, sd_hetero, rho_base, mu_pop,
                              "lavaan only v4 condshrink power only",
                              "cells C & D")
  mlm_res_alt_AB <- run_mlm_pass(file.path(Data, sprintf("mlm results %s.rds",
                                                         mlm_tag_alt_AB)),
                                 rep_fn         = run_one_rep_mlm_power_only_AB,
                                 checkpoint_tag = sprintf("mlm alt %s", 
                                                          mlm_tag_alt_AB),
                                 label          = "appendix MLM power-only (cells A & B)",
                                 force_fresh    = force_rerun_MLM,
                                 dist_A_neg_arg = dist_A_true_ide_alt,
                                 dist_B_neg_arg = dist_B_true_ide_alt)
  mlm_res_alt_CD <- run_mlm_pass(file.path(Data, sprintf("mlm results %s.rds", 
                                                         mlm_tag_alt_CD)),
                                   rep_fn          = run_one_rep_mlm_power_only_CD,
                                   checkpoint_tag  = sprintf("mlm alt %s", 
                                                             mlm_tag_alt_CD),
                                   label           = "appendix MLM power-only (cells C & D)",
                                   force_fresh     = force_rerun_MLM,
                                   pool_set        = appendix_pool_set)
  mlm_res_alt    <- merge(mlm_res_alt_AB, mlm_res_alt_CD, by = "rep")
  log_step("Appendix module verified, cells A to D, power only (alternative tracking-preserving true_ide; Type 1 error reused from primary Table 2)") }

## PROGRAM    -- TARGET THRESHOLD SWEEP (K REPLICATIONS, CELLS C & D) ----------

log_step("TARGET THRESHOLD SWEEP (K REPLICATIONS, CELLS C & D)", 1)

# Extends primary analysis (z = 1.645, top 5%) across 4 threshold levels evenly
# spaced on the log scale of tail probability (25%, 5%, 1%, 0.2%). Only relevant
# for C & D (as both involve targeted allocation on baseline); Cells A & B are
# invariant to threshold. Each threshold level has its own pool set (six pools:
# control C / D + treat C / D x no_ide / true_ide).

# sweep_K: this sweep (both the Oldham & MLM legs below) runs at
# min(K, sweep_K_max) (PART ONE), so the main Table 2 passes & the PART FIVE
# diagnostic keep the full K. Both legs - especially the MLM leg - are the
# most expensive part of the whole script (each threshold can run for many
# hours at large K). sweep_tag/sweep_mlm_tag below & the severity plot file
# name all carry sweep_K.

sweep_K        <- min(K, sweep_K_max)
sweep_tag      <- cache_tag(cache_version, sweep_K, sweep_pct, mu_pop, sd_homo,
                            sd_hetero, rho_base, mean_mult, sd_mult_true_ide, Nsmp)
sweep_res_file <- file.path(Sweep_Data, sprintf("sweep results %s.rds", sweep_tag))

# One pool set (six truncated-distribution specifications, as in BUILD
# TARGETED-ALLOCATION POOLS above) per threshold, held in sweep_pools for the
# Oldham & MLM legs below & PART FIVE's diagnostic

sweep_pools <- list()
for (pct_i in seq_along(sweep_pct)) {
  pct   <- sweep_pct[pct_i]
  z     <- sweep_z[pct_i]
  z_tag <- sprintf("z%04d", round(z * 1000))
  ps    <- list()
  for (spec_row in seq_len(nrow(pool_spec))) {
    spec      <- pool_spec[spec_row, ]
    sd_vec    <- if (spec$sd_regime == "hetero") sd_hetero else sd_homo
    pool_name <- if (spec$role == "control") 
      sprintf("pool_control_%s", spec$cell) else
      sprintf("pool_treat_%s_%s", spec$cell, spec$ide_type)
    if (spec$role == "control") {
      pool    <- build_targeted_pool(mu_vec = mu_pop, sd_vec = sd_vec,
                                     rho_base = rho_base, ref_mean = mu_pop, 
                                     ref_sd = sd_vec, z_cut = z)
    } else {
      dist    <- build_treatment_distribution(mu_vec = mu_pop, sd_vec = sd_vec, 
                                              rho_base = rho_base, 
                                              mean_mult = mean_mult, 
                                              sd_mult_true_ide = sd_mult_true_ide,
                                              ide_type = spec$ide_type)
      pool    <- build_targeted_pool(mu_vec = dist$mu_vec, sd_vec = sd_vec,
                                     rho_base = rho_base, ref_mean = mu_pop, 
                                     ref_sd = sd_vec, 
                                     Sigma_override = dist$Sigma, z_cut = z) }
    ps[[pool_name]] <- pool }
  sweep_pools[[z_tag]] <- list(
    pct           = pct,
    z             = z,
    z_tag         = z_tag,
    pool_set      = list(
      treat_C     = ps[["pool_treat_C_no_ide"]],
      treat_C_neg = ps[["pool_treat_C_true_ide"]],
      control_C   = ps[["pool_control_C"]],
      treat_D     = ps[["pool_treat_D_no_ide"]],
      treat_D_neg = ps[["pool_treat_D_true_ide"]],
      control_D   = ps[["pool_control_D"]])) }

# PER-THRESHOLD CACHING: each threshold's summary row is cached to its own
# file (keyed on sweep_tag + that threshold's z_tag) & written to disk the
# moment that threshold finishes, rather than only once at the very end of
# the full 4-threshold loop. Without this, killing the script partway
# through threshold 3 of 4 would lose thresholds 1 & 2 as well - their
# per-replication checkpoint files self-delete on completion (see
# run_reps_chunked()), so the ONLY durable record of a finished
# threshold used to be the final sweep_summary object, which was never
# written until every threshold was done. Now a completed threshold survives
# a kill/crash on its own & is skipped (not recomputed) on the next run.

log_step(sprintf("TARGET THRESHOLD SWEEP: K = %s REPLICATIONS PER THRESHOLD",
                 format(as.integer(sweep_K), big.mark = ",")))
sweep_results  <- vector("list", length(sweep_pools))
for (pool_i in seq_along(sweep_pools)) {
  sp             <- sweep_pools[[pool_i]]
  sweep_row_file <- file.path(Sweep_Data, sprintf("sweep results %s %s.rds", 
                                                  sweep_tag, sp$z_tag))
  if (!force_rerun_SWEEP && file.exists(sweep_row_file)) {
    sweep_results[[pool_i]] <- readRDS(sweep_row_file)
    log_step(sprintf("Threshold %g%% (z = %.3f) -- loaded from cache", 
                     sp$pct, sp$z))
  } else {
    log_step(sprintf("Threshold %g%% (z = %.3f) ...", sp$pct, sp$z))
    t_sw <- Sys.time()
    raw  <- run_reps_chunked(sweep_K, parallel_K,
                             chunk_size     = chunk_size,
                             rep_fn         = run_one_replication_CD,
                             cl             = get_script_cl(),
                             checkpoint_tag = sprintf("sweep %s %s", sweep_tag, 
                                                      sp$z_tag),
                             force_fresh    = force_rerun_SWEEP,
                             data_dir       = Sweep_Data,
                             z_tag          = sp$z_tag)
    log_done(t_sw, sprintf("%g%% threshold, K = %s", sp$pct, 
                           format(as.integer(sweep_K), big.mark = ",")))
    sweep_row <- data.frame(
      threshold_pct      = sp$pct,
      z_cut              = sp$z,
      Old_C_T1e_NoComp   = reject_rate(raw$p_C_int_only),
      Old_C_T1e_Comp     = reject_rate(raw$p_C_vs_comparator),
      Old_C_Power_NoComp = power_rate(raw$p_C_int_only_neg),
      Old_C_Power_Comp   = power_rate(raw$p_C_vs_comparator_neg),
      Old_D_T1e_NoComp   = reject_rate(raw$p_D_int_only),
      Old_D_T1e_Comp     = reject_rate(raw$p_D_vs_comparator),
      Old_D_Power_NoComp = power_rate(raw$p_D_int_only_neg),
      Old_D_Power_Comp   = power_rate(raw$p_D_vs_comparator_neg))
    saveRDS(sweep_row, sweep_row_file)
    sweep_results[[pool_i]] <- sweep_row } }
sweep_summary <- do.call(rbind, sweep_results)
saveRDS(sweep_summary, sweep_res_file)
log_step(sprintf("Target threshold sweep (Oldham leg) complete -- results cached to %s", 
                 sweep_res_file))
log_step("Target threshold sweep summary (Oldham leg):"); print(sweep_summary)

## PROGRAM    -- TARGET THRESHOLD SWEEP, MLM LEG (CELLS C & D) ----------------

# Runs the SAME four thresholds through the MLM (lavaan) route, on the
# identical draws as the Oldham leg above (run_one_replication_CD_mlm() uses
# the same seed_for() keys as run_one_replication_CD()), so the two methods
# are directly comparable point-for-point across the severity sweep. Shares
# sweep_K with the Oldham leg above, min(K, sweep_K_max) - set K low (e.g.
# 100) to smoke-test the whole pipeline quickly, then raise it for a real run,
# since the "MLM stays near 5%" claim only firms up at a larger sweep_K.

sweep_mlm_tag      <- cache_tag(cache_version, sweep_K, sweep_pct, mu_pop, 
                                sd_homo, sd_hetero, rho_base, mean_mult, 
                                sd_mult_true_ide, Nsmp, "mlm leg", "conv counts")
sweep_mlm_res_file <- file.path(Sweep_Data, sprintf("sweep mlm results %s.rds", 
                                                    sweep_mlm_tag))

# PER-THRESHOLD CACHING (see matching note on the Oldham leg above): each
# threshold's row is cached & written to disk as soon as that threshold
# finishes, not only at the very end of the 4-threshold loop - this leg is
# by far the most expensive part of the whole script (each threshold can run
# several hours), so losing a completed threshold to a kill/crash because
# only the FINAL aggregate was ever persisted was the single biggest risk in
# the entire pipeline. A completed threshold's row file now survives on its
# own & is skipped (not recomputed) on the next run. Each threshold's MLM
# exclusion counts (mlm_exclusions(), PART TWO) are cached beside its row file
# in the same way, since the raw per-replication frame is not kept.

log_step(sprintf("TARGET THRESHOLD SWEEP (MLM LEG): K = %s REPLICATIONS PER THRESHOLD",
                 format(as.integer(sweep_K), big.mark = ",")))
sweep_mlm_results  <- vector("list", length(sweep_pools))
sweep_mlm_excl     <- vector("list", length(sweep_pools))
for (pool_i in seq_along(sweep_pools)) {
  sp                  <- sweep_pools[[pool_i]]
  sweep_mlm_row_file  <- file.path(Sweep_Data, 
                                   sprintf("sweep mlm results %s %s.rds", 
                                           sweep_mlm_tag, sp$z_tag))
  sweep_mlm_excl_file <- file.path(Sweep_Data, 
                                   sprintf("sweep mlm exclusions %s %s.rds", 
                                           sweep_mlm_tag, sp$z_tag))
  if (!force_rerun_SWEEP && file.exists(sweep_mlm_row_file) && 
      file.exists(sweep_mlm_excl_file)) {
    sweep_mlm_results[[pool_i]] <- readRDS(sweep_mlm_row_file)
    sweep_mlm_excl[[pool_i]]    <- readRDS(sweep_mlm_excl_file)
    log_step(sprintf("MLM threshold %g%% (z = %.3f) -- loaded from cache", sp$pct, sp$z))
  } else {
    log_step(sprintf("MLM threshold %g%% (z = %.3f) ...", sp$pct, sp$z))
    t_sw <- Sys.time()
    raw  <- run_reps_chunked(sweep_K, parallel_K = 1, 
                             chunk_size     = mlm_chnk_size,
                             rep_fn         = run_one_replication_CD_mlm,
                             cl             = get_script_cl(needs_parallel = TRUE),
                             checkpoint_tag = sprintf("sweep mlm %s %s", 
                                                      sweep_mlm_tag, sp$z_tag),
                             force_fresh    = force_rerun_SWEEP,
                             data_dir       = Sweep_Data,
                             z_tag          = sp$z_tag)
    log_done(t_sw, sprintf("MLM %g%% threshold, K = %s", sp$pct, 
                           format(as.integer(sweep_K), big.mark = ",")))
    sweep_mlm_row <- data.frame(
      threshold_pct      = sp$pct,
      z_cut              = sp$z,
      MLM_C_T1e_NoComp   = reject_rate(raw$p_C_int_only),
      MLM_C_T1e_Comp     = reject_rate(raw$p_C_vs_comparator),
      MLM_C_Power_NoComp = power_rate(raw$p_C_int_only_neg),
      MLM_C_Power_Comp   = power_rate(raw$p_C_vs_comparator_neg),
      MLM_D_T1e_NoComp   = reject_rate(raw$p_D_int_only),
      MLM_D_T1e_Comp     = reject_rate(raw$p_D_vs_comparator),
      MLM_D_Power_NoComp = power_rate(raw$p_D_int_only_neg),
      MLM_D_Power_Comp   = power_rate(raw$p_D_vs_comparator_neg))
    sweep_mlm_excl[[pool_i]] <- mlm_exclusions(raw, sprintf("Sweep %g%%", sp$pct))
    saveRDS(sweep_mlm_excl[[pool_i]], sweep_mlm_excl_file)
    saveRDS(sweep_mlm_row, sweep_mlm_row_file)
    sweep_mlm_results[[pool_i]] <- sweep_mlm_row } }
sweep_mlm_summary <- do.call(rbind, sweep_mlm_results)
saveRDS(sweep_mlm_summary, sweep_mlm_res_file)
log_step(sprintf("Target threshold sweep (MLM leg) complete -- results cached to %s", 
                 sweep_mlm_res_file))
log_step("Target threshold sweep summary (MLM leg):"); print(sweep_mlm_summary)

# Merge Oldham & MLM legs into one summary (by threshold_pct/z_cut, which are
# identical across both legs since both are built from sweep_pct/sweep_z) so
# the severity plot below reads from a single object

sweep_summary <- merge(sweep_summary, sweep_mlm_summary, 
                       by = c("threshold_pct", "z_cut"))
sweep_summary <- sweep_summary[order(-sweep_summary$threshold_pct), ]
rownames(sweep_summary) <- NULL

## PROGRAM    -- PLOT (OLDHAM VS MLM, WITH-COMPARATOR TYPE 1 ERROR & POWER) ----

# Shows the paper's core threshold-severity argument in one figure: as the
# targeted-recruitment cut moves from 25% (mild truncation) to 0.2% (severe
# truncation), Oldham's Fisher z contrast is increasingly stretched beyond
# the assumptions its reference variance relies on (see PART FIVE's
# diagnostic: the empirical/nominal variance ratio grows with severity), so
# its "with comparator" Type 1 error climbs well above the nominal 5% line.
# The MLM route stays near 5% because the growth-model LRT does not lean on
# that same fixed reference-variance formula. Only the "with comparator" test
# is plotted (not "intervention-only", which sits pinned at 100% for both
# methods under targeted recruitment regardless of severity -- see Table 2 --
# and would swamp the y-axis).
# Type 1 error & power are combined into ONE chart (facet_grid: metric as
# rows - Type 1 error above Power - cell as columns), rather than two
# separate figures, so severity's effect on both is visible together in one
# place. sweep_K is annotated in the subtitle since
# the MLM leg's proximity to 5% is a claim that only firms up as sweep_K
# grows (small sweep_K leaves MLM's own Monte Carlo noise visible around the
# 5% line) -- see the sweep_K note above (PART FOUR), min(K, sweep_K_max).
# X-AXIS: log-scale tail probability, REVERSED so severity increases
# left-to-right (25% mild truncation on the left, 0.2% severe truncation on
# the right) - the more intuitive reading direction for a "severity" axis.
# Implemented by plotting log10(threshold_pct) directly & using
# scale_x_reverse() with breaks/labels supplied in that same log10 space,
# rather than relying on trans = c("log10","reverse") composition (version-
# dependent across ggplot2/scales releases) - this way is simple & robust.

severity_power_data <- rbind(
  data.frame(threshold_pct = sweep_summary$threshold_pct, cell = "C", 
             method = "Oldham", metric = "Type 1 error", 
             value = sweep_summary$Old_C_T1e_Comp),
  data.frame(threshold_pct = sweep_summary$threshold_pct, cell = "D", 
             method = "Oldham", metric = "Type 1 error", 
             value = sweep_summary$Old_D_T1e_Comp),
  data.frame(threshold_pct = sweep_summary$threshold_pct, cell = "C", 
             method = "MLM", metric = "Type 1 error", 
             value = sweep_summary$MLM_C_T1e_Comp),
  data.frame(threshold_pct = sweep_summary$threshold_pct, cell = "D", 
             method = "MLM", metric = "Type 1 error", 
             value = sweep_summary$MLM_D_T1e_Comp),
  data.frame(threshold_pct = sweep_summary$threshold_pct, cell = "C", 
             method = "Oldham", metric = "Power", 
             value = sweep_summary$Old_C_Power_Comp),
  data.frame(threshold_pct = sweep_summary$threshold_pct, cell = "D", 
             method = "Oldham", metric = "Power", 
             value = sweep_summary$Old_D_Power_Comp),
  data.frame(threshold_pct = sweep_summary$threshold_pct, cell = "C",
             method = "MLM", metric = "Power", 
             value = sweep_summary$MLM_C_Power_Comp),
  data.frame(threshold_pct = sweep_summary$threshold_pct, cell = "D",
             method = "MLM", metric = "Power", 
             value = sweep_summary$MLM_D_Power_Comp))
severity_power_data <- transform(severity_power_data, 
                                 cell   = factor(cell, levels = c("C", "D"),
                                                 labels = c("Scenario C (constant variance)",
                                                            "Scenario D (non-constant variance)")),
                                 method = factor(method, 
                                                 levels = c("Oldham", "MLM")),
                                 metric = factor(metric, 
                                                 levels = c("Type 1 error",
                                                            "Power")),
                                 log_threshold = log10(threshold_pct))

# 5% reference line only on the Type 1 error row (power is not expected to
# sit near 5%) - a small per-facet data frame so geom_hline() only draws in
# the matching facet panels under facet_grid()

hline_data <- data.frame(metric = factor("Type 1 error", 
                                         levels = c("Type 1 error", "Power")),
                         yintercept = 5)
severity_power_plot <- ggplot(severity_power_data, aes(x = log_threshold, 
                                                       y = value, 
                                                       colour = method)) +
  geom_hline(data = hline_data, aes(yintercept = yintercept), 
             linetype = "dashed", colour = "grey50") +
  geom_line() + geom_point(size = 2) +
  scale_x_reverse(breaks = log10(sweep_pct), 
                  labels = sprintf("%g%%", sweep_pct)) +
  scale_colour_manual(values = c(Oldham = "#D55E00", MLM = "#0072B2")) +
  facet_grid(metric ~ cell, scales = "free_y") +
  # No title/subtitle: this plot's title belongs in the Word manuscript, not
  # baked into the image - see the note on SaveJpg() below re: image height.
  labs(x      = "Recruitment threshold (top tail %, log scale, mild to severe)",
       y      = "Rate (%)",
       colour = "Method") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

# A4 width (8.27 in) by 40% of A4 height (4.68 in), the shape of Figure S1 as
# placed in the Word manuscript (29 July script); the 21 Sept copy had reverted
# to 90% of A4 height, a tall portrait image that no longer matched Word

SaveJpg(severity_power_plot, 
        sprintf("Sweep_type1 error & power vs severity, Oldham vs MLM sweep_K%d", 
                sweep_K), width = 8.27, height = 11.69 * 0.4)
log_step("Severity plot (Oldham vs. MLM, with-comparator Type 1 error & power, stacked) saved in Plots")

## #############################################################################
## PART FIVE  -- CONFIRMATORY DIAGNOSTIC (FISHER Z REFERENCE VARIANCE)
## PROGRAM    -- POOLS FOR DIAGNOSTICS -----------------------------------------

log_step("CONFIRMATORY DIAGNOSTIC: POOLS", 1)

# Takes the no_ide treatment & control pools for the 5% & 0.2% thresholds
# from the sweep's own pool sets (sweep_pools, PART FOUR); diag_z_tags fixes
# which two thresholds are diagnosed

diag_z_tags    <- c("z1645", "z2878")
diag_pools     <- list()
for (z_tag in diag_z_tags) {
  ps           <- sweep_pools[[z_tag]]$pool_set
  diag_pools[[z_tag]] <- list(
    treat_C    = ps$treat_C,
    control_C  = ps$control_C,
    treat_D    = ps$treat_D,
    control_D  = ps$control_D) }

# Runs K replications at EACH of the two thresholds in diag_z_tags, using
# run_one_replication_CD_z() (PART TWO) rather than run_one_replication_CD():
# the sweep only needed a p-value per replication, but this diagnostic needs
# the raw per-arm Fisher z values themselves (z1_C, z0_C, z1_D, z0_D, plus
# each arm's n) so diagnose_cell() below can compare their EMPIRICAL variance
# against the NOMINAL variance formula the z-test assumes. Results across
# both thresholds are stacked into one data frame (z_tag column identifies
# which threshold each row came from) & cached as a single object, since
# both thresholds share the same K & the same downstream summary step.

diag_tag       <- cache_tag(cache_version, K, diag_z_tags, sweep_tag)
diag_res_file  <- file.path(Sweep_Data, 
                            sprintf("diagnostic results %s.rds", diag_tag))
if (force_rerun_SWEEP || !file.exists(diag_res_file)) {
  diag_results <- vector("list", length(diag_z_tags))
  for (tag_i in seq_along(diag_z_tags)) {
    z_tag      <- diag_z_tags[tag_i]
    t_dg       <- Sys.time()
    log_step(sprintf("Diagnostic replications for threshold %s ...", z_tag))
    raw        <- run_reps_chunked(K, parallel_K,
                                   chunk_size     = chunk_size,
                                   rep_fn         = run_one_replication_CD_z,
                                   cl             = get_script_cl(),
                                   checkpoint_tag = sprintf("diag %s %s", 
                                                            diag_tag, z_tag),
                                   force_fresh    = force_rerun_SWEEP,
                                   data_dir       = Sweep_Data,
                                   z_tag          = z_tag)
    log_done(t_dg, 
             sprintf("threshold %s, K = %s", z_tag, format(as.integer(K),
                                                           big.mark = ",")))
    raw$z_tag <- z_tag
    diag_results[[tag_i]] <- raw }
  diag_raw    <- do.call(rbind, diag_results)
  saveRDS(diag_raw, diag_res_file)
  log_step(sprintf("Diagnostic replications complete -- cached to %s", 
                   diag_res_file)) } else {
    diag_raw  <- readRDS(diag_res_file)
    log_step("Diagnostic replications loaded from cache") }

# Loops over both thresholds (diag_z_tags) x both cells (C, D) - 4 combinations
# - calling diagnose_cell() (PART TWO) once per combination on that
# combination's slice of diag_raw. Collects a one-row summary per combination
# (diag_summary_rows, feeding Table A2 & the printed table below), row-bound
# into a single data frame once the loop finishes.

diag_summary_rows <- list()
for (z_tag in diag_z_tags) {
  sub_raw         <- diag_raw[diag_raw$z_tag == z_tag, ]
  dC              <- diagnose_cell(sub_raw$z1_C, sub_raw$z0_C, sub_raw$n1_C, 
                                   sub_raw$n0_C)
  dD              <- diagnose_cell(sub_raw$z1_D, sub_raw$z0_D, sub_raw$n1_D, 
                                   sub_raw$n0_D)
  diag_summary_rows[[paste0(z_tag, "_C")]] <- data.frame(
    z_tag = z_tag, cell = "C", mean_contrast = dC$mean_contrast, 
    empirical_variance = dC$empirical_variance, 
    nominal_variance = dC$nominal_variance, ratio = dC$ratio)
  diag_summary_rows[[paste0(z_tag, "_D")]] <- data.frame(
    z_tag = z_tag, cell = "D", mean_contrast = dD$mean_contrast,
    empirical_variance = dD$empirical_variance, 
    nominal_variance = dD$nominal_variance, ratio = dD$ratio) }
diag_summary           <- do.call(rbind, diag_summary_rows)
rownames(diag_summary) <- NULL
print(diag_summary)
log_step("Interpretation guide: mean_contrast near zero indicates no location bias")
log_step("(arms agree in expectation); ratio above 1 indicates the true variance")
log_step("of the Fisher z contrast exceeds the nominal formula, which is exactly")
log_step("the miscalibration that would inflate Type 1 error in the z-test")

## PROGRAM    -- STOP SCRIPT-WIDE CLUSTER (NOTHING PARALLEL RUNS AFTER THIS) --

if (!is.null(script_cl_handle) ) {
  stopCluster(script_cl_handle)
  log_step("Script-wide cluster stopped") 
} else if (script_cl_built) {
  log_step("No script-wide cluster was needed (K <= parallel_K at the top level)") 
} else { log_step("No script-wide cluster was ever built (every stage loaded from cache)") }

## #############################################################################
## PART SIX   -- FIGURES & TABLES (TYPE 1 ERROR / POWER, NO_IDE / true_ide)
## PROGRAM    -- TABLE 2 (TYPE 1 ERROR / POWER, OLDHAM / MLM) ------------------

# Table 2, long format: one row per cell x IDE status (None = no_ide, -ve =
# true_ide), assembled via table2_row() (PART TWO, PLOT SAVE & TABLE FORMAT
# HELPERS). Every row carries all four value columns (Oldham T1 Error, Oldham
# Power, MLM T1 Error, MLM Power), but only the column matching that row's
# IDE status is populated - the others show "-" (not NA/blank), since a Type
# 1 error is only defined on a no_ide draw & a power figure is only defined
# on a true_ide draw; showing "-" makes that explicit rather than leaving an
# ambiguous empty cell. Intervention Arm block: cells A-D, Oldham's 
# intervention-only test vs the MLM single-arm analogue (lavaan_one_arm(), 
# PART TWO) - both routes now populated, evaluated on identical draws (PART 
# FOUR, run_one_replication_mlm_AB() / _CD())

intervention_rows <- rbind(
  table2_row("Intervention Arm", "A", "None", TRUE,
             oldham_t1 = reject_rate(verify_res$p_A_int_only), 
             mlm_t1 = reject_rate(mlm_res$p_A_int_only) ),
  table2_row("Intervention Arm", "A", "-ve",  FALSE,
             oldham_pow = power_rate(verify_res$p_A_int_only_neg), 
             mlm_pow = power_rate(mlm_res$p_A_int_only_neg) ),
  table2_row("Intervention Arm", "B", "None", TRUE,
             oldham_t1 = reject_rate(verify_res$p_B_int_only), 
             mlm_t1 = reject_rate(mlm_res$p_B_int_only) ),
  table2_row("Intervention Arm", "B", "-ve",  FALSE,
             oldham_pow = power_rate(verify_res$p_B_int_only_neg), 
             mlm_pow = power_rate(mlm_res$p_B_int_only_neg) ),
  table2_row("Intervention Arm", "C", "None", TRUE,
             oldham_t1 = reject_rate(verify_res$p_C_int_only), 
             mlm_t1 = reject_rate(mlm_res$p_C_int_only) ),
  table2_row("Intervention Arm", "C", "-ve",  FALSE,
             oldham_pow = power_rate(verify_res$p_C_int_only_neg), 
             mlm_pow = power_rate(mlm_res$p_C_int_only_neg) ),
  table2_row("Intervention Arm", "D", "None", TRUE,
             oldham_t1 = reject_rate(verify_res$p_D_int_only), 
             mlm_t1 = reject_rate(mlm_res$p_D_int_only) ),
  table2_row("Intervention Arm", "D", "-ve",  FALSE,
             oldham_pow = power_rate(verify_res$p_D_int_only_neg), 
             mlm_pow = power_rate(mlm_res$p_D_int_only_neg) ) )

# With Comparator block: cells B-D only (Cell A has no comparator arm);
# Oldham & MLM both populated, evaluated on the same draws (PART FOUR)

comparator_rows <- rbind(
  table2_row("With Comparator", "B", "None", TRUE,
             oldham_t1 = reject_rate(verify_res$p_B_vs_comparator), 
             mlm_t1 = reject_rate(mlm_res$p_B_vs_comparator) ),
  table2_row("With Comparator", "B", "-ve",  FALSE,
             oldham_pow = power_rate(verify_res$p_B_vs_comparator_neg), 
             mlm_pow = power_rate(mlm_res$p_B_vs_comparator_neg) ),
  table2_row("With Comparator", "C", "None", TRUE,
             oldham_t1 = reject_rate(verify_res$p_C_vs_comparator), 
             mlm_t1 = reject_rate(mlm_res$p_C_vs_comparator) ),
  table2_row("With Comparator", "C", "-ve",  FALSE,
             oldham_pow = power_rate(verify_res$p_C_vs_comparator_neg), 
             mlm_pow = power_rate(mlm_res$p_C_vs_comparator_neg) ),
  table2_row("With Comparator", "D", "None", TRUE,
             oldham_t1 = reject_rate(verify_res$p_D_vs_comparator), 
             mlm_t1 = reject_rate(mlm_res$p_D_vs_comparator) ),
  table2_row("With Comparator", "D", "-ve",  FALSE,
             oldham_pow = power_rate(verify_res$p_D_vs_comparator_neg), 
             mlm_pow = power_rate(mlm_res$p_D_vs_comparator_neg) ) )

Table2 <- rbind(intervention_rows, comparator_rows)
log_step("TABLE 2: ASSEMBLED", 1); print(Table2)

## PROGRAM    -- TABLE A1 (ALTERNATIVE TRACKING-PRESERVING true_ide, POWER ONLY) --

# Power-only appendix table (tableA1_row(), PART TWO), built from
# verify_res_alt/mlm_res_alt (PART FOUR) - the alternative true_ide definition
# (build_treat_dist_condshrink(), PART TWO) kept as an Appendix
# robustness check against the primary tracking-weakening assumption. Type 1
# error is NOT repeated here (identical to Table 2 under both true_ide
# specifications, since no_ide is unaffected by which true_ide builder is
# used), & the intervention-arm-only block is restricted to Cell A only -
# unlike Table2, Cells B-D are omitted from Intervention Arm Only entirely,
# since that evaluation is invalid there regardless of which true_ide
# specification is in play, matching the paper's own Table A1. If
# skip_appendix is TRUE, verify_res_alt/mlm_res_alt were never computed (PART
# FOUR), so TableA1 is set to NULL instead - OUTPUT ALL RESULTS below skips
# writing its xlsx sheet in that case.

if (skip_appendix) {
  TableA1 <- NULL
  log_step("skip_appendix = TRUE -- Table A1 omitted", 1)
} else {
  appendix_intervention_rows <- rbind(
    tableA1_row("Intervention Arm", "A", "None", TRUE),
    tableA1_row("Intervention Arm", "A", "-ve",  FALSE,
                oldham_pow = power_rate(verify_res_alt$p_A_int_only_neg), 
                mlm_pow = power_rate(mlm_res_alt$p_A_int_only_neg) ) )
  appendix_comparator_rows <- rbind(
    tableA1_row("With Comparator", "B", "None", TRUE),
    tableA1_row("With Comparator", "B", "-ve",  FALSE,
                oldham_pow = power_rate(verify_res_alt$p_B_vs_comparator_neg), 
                mlm_pow = power_rate(mlm_res_alt$p_B_vs_comparator_neg) ),
    tableA1_row("With Comparator", "C", "None", TRUE),
    tableA1_row("With Comparator", "C", "-ve",  FALSE,
                oldham_pow = power_rate(verify_res_alt$p_C_vs_comparator_neg), 
                mlm_pow = power_rate(mlm_res_alt$p_C_vs_comparator_neg) ),
    tableA1_row("With Comparator", "D", "None", TRUE),
    tableA1_row("With Comparator", "D", "-ve",  FALSE,
                oldham_pow = power_rate(verify_res_alt$p_D_vs_comparator_neg), 
                mlm_pow = power_rate(mlm_res_alt$p_D_vs_comparator_neg) ) )
  TableA1 <- rbind(appendix_intervention_rows, appendix_comparator_rows)
  log_step("TABLE A1: ASSEMBLED (alternative tracking-preserving true_ide, power only)",
           1); print(TableA1) }

## PROGRAM    -- MLM EXCLUSIONS (NO USABLE P-VALUE, NOT CONVERGED) -------------

# reject_rate() & power_rate() drop NA p-values, so every MLM rate in Table 2,
# Table A1 & the severity sweep is computed over fewer than K replications
# wherever a fit returned no usable p-value. MLMExclusions stacks the
# mlm_exclusions() counts (PART TWO) for the Table 2 pass (mlm_res), the Table
# A1 pass (mlm_res_alt, unless skip_appendix) & every sweep threshold (MLM
# leg, PART FOUR); the maximum no usable p-value count is what the manuscript
# quotes as excluded replications

MLMExclusions <- rbind(
  mlm_exclusions(mlm_res, "Table 2"),
  if (!skip_appendix) mlm_exclusions(mlm_res_alt, "Table A1"),
  do.call(rbind, sweep_mlm_excl))
log_step("MLM EXCLUSIONS: REPLICATIONS WITH NO USABLE P-VALUE / NOT CONVERGED", 1)
print(MLMExclusions)
log_step(sprintf("Maximum no usable p-value count: %d   Maximum not converged count: %d",
                 max(MLMExclusions$`No usable p-value`), 
                 max(MLMExclusions$`Not converged`)))

## PROGRAM    -- OUTPUT ALL RESULTS --------------------------------------------

# Table A2 (paper format): diag_summary (PART FIVE) holds raw-scale variance
# figures & the raw z_tag/cell keys; the paper reports Table A2 with Threshold/
# Scenario labels & every variance figure x10-3, so that relabelling &
# rescaling happens here, right at output, rather than upstream where
# diag_summary is also used unscaled by the interpretation-guide log_steps.

threshold_label <- c(z1645 = "5%", z2878 = "0.2%")
scenario_label  <- c(C = "C, constant variance", D = "D, non-constant variance")
TableA2 <- data.frame(
  Threshold                    = threshold_label[diag_summary$z_tag],
  Scenario                     = scenario_label[diag_summary$cell],
  `Mean Contrast (x10-3)`      = round(diag_summary$mean_contrast*1000, 2),
  `Reference Variance (x10-3)` = round(diag_summary$nominal_variance*1000, 2),
  `Empirical Variance (x10-3)` = round(diag_summary$empirical_variance*1000, 2),
  Ratio                        = round(diag_summary$ratio, 2),
  check.names = FALSE, stringsAsFactors = FALSE)

# xlsx file is written to Data/K_<K>/ (alongside that K's cached .rds
# artefacts) under a fixed name - NOT dated (neither by Sys.Date() nor by the
# script file's own name/date) - since the K-specific subfolder & its
# unique-tag-named .rds files already identify a run; a date in the filename
# on top of that added no information & only caused confusion between which
# run's parameters actually produced a given dated file.

xlsx_path <- file.path(Data, "IDE results tables.xlsx")
wb        <- createWorkbook()
addWorksheet(wb, "Table 2")
addWorksheet(wb, "Table A2")
writeData(wb, "Table 2", Table2)
writeData(wb, "Table A2", TableA2)

# Table A1 sheet omitted entirely when skip_appendix = TRUE (PART ONE) rather
# than added with an empty/NULL TableA1, so a skipped run's workbook doesn't
# carry a misleadingly-present but blank sheet

if (!skip_appendix) {
  addWorksheet(wb, "Table A1")
  writeData(wb, "Table A1", TableA1) }

# Figure A1 sheet: the merged sweep_summary (PART FOUR) behind the severity
# plot, both legs, every threshold; MLM exclusions sheet: MLMExclusions above

addWorksheet(wb, "Figure A1")
writeData(wb, "Figure A1", sweep_summary)
addWorksheet(wb, "MLM exclusions")
writeData(wb, "MLM exclusions", MLMExclusions)
saveWorkbook(wb, xlsx_path, overwrite = TRUE)
log_step(sprintf("Results workbook written: %s", xlsx_path) )

## #############################################################################
## END