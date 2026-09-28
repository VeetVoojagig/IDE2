IDE2 – EVALUATING CHANGE IN HETEROSCEDASTIC DATA
================================================================================

Author: Mark S Gilthorpe


WHAT THIS PROJECT IS ABOUT
--------------------------------------------------------------------------------

This project is a Monte Carlo simulation study, written in R, that looks at a
common question in evaluation research: did an intervention change how much
people differ from one another over time, as well as changing their average?

This is called an intervention differential effect (IDE). A negative IDE means
the intervention weakens "tracking", i.e. how strongly a child's early weight
predicts their later weight. The usual tool for testing this is Oldham's
method. It correlates the change between baseline and follow-up with the
average of the two, and asks whether that correlation differs from zero.

Oldham's method assumes the outcome's variance stays the same over time, and
that allocation to the intervention does not depend on the baseline value.
Real data often break both assumptions:

  1. Variance often grows with time (heteroscedasticity). Children's weight
     spreads out as they get older, for example.
  2. Interventions are often targeted at people with high baseline values,
     e.g. only children above an "at risk" weight are recruited. This
     truncates the baseline distribution.

The project shows how each of these problems, alone and together, makes
standard tests of change misleading. It also shows what a valid comparison
has to look like.


THE 2 x 2 SIMULATION FRAMEWORK
--------------------------------------------------------------------------------

The study crosses two factors:

```
  allocation   : random, or targeted on baseline (truncation)
  variance     : homoscedastic (constant SD 1.0 / 1.0 / 1.0), or
                 heteroscedastic (rising SD 1.0 / 1.5 / 2.0)
```

This gives four cells, called Scenarios A to D in the manuscript:

```
  A  homoscedastic,   random    Positive control / calibration check.
                                Intervention arm only; no comparator.
  B  heteroscedastic, random    A population-wide measure, e.g. a sugar levy
                                or advertising ban. The paper's central
                                empirical result.
  C  homoscedastic,   targeted  An extreme version of Beggs et al.; shows the
                                effect of very severe truncation.
  D  heteroscedastic, targeted  Both problems together, e.g. a letter
                                sent to children above a weight threshold.
```

The data are girls' weight at ages 1, 2.5 and 4 years. Means come from
CDC/LMS weight-for-age reference data (9.67, 12.78 and 15.88 kg). The lag
correlation between occasions is 0.9 raised to the time gap in years. Each
arm has N = 500.

Every cell is run under two conditions:

```
  no_ide   The intervention shifts the mean only (Type 1 error is measured
           here).
  true_ide The intervention also weakens tracking by shrinking the follow-up
           SD, reaching 0.875 of the control SD by age 4 (power is measured
           here).
```

The intervention effect is built into the distribution that the treatment arm
is drawn from. It is never added to the data afterwards.

For the targeted cells (C and D), each arm is drawn afresh in every
replication from its own baseline-truncated distribution. Both arms are cut at
the same baseline threshold (top 5% by default, z > 1.645). The comparator has
to come from the same truncated population as the treatment arm, otherwise the
two groups would not be comparable.


METHODS COMPARED
--------------------------------------------------------------------------------

Each simulated dataset is analysed in two ways:

```
  Oldham   Oldham's method, as described above.
  MLM      A latent growth model fitted in lavaan, with intercept and slope
           centred at age 2.5. It tests whether the intercept-slope
           covariance differs between arms, using a likelihood ratio test.
```

Each method is also evaluated in two ways:

```
  Intervention arm only    The intervention group alone, tested against no
                           change.
  With comparator          Intervention vs. control group.
```


WHAT THE PROJECT ACHIEVES
--------------------------------------------------------------------------------

  - It measures Type 1 error and power for every combination of cell, method
    and evaluation, over K = 1,000,000 replications.
  - It shows that testing the intervention arm alone breaks down whenever
    variance is not constant (Cells B to D), even when there is no real IDE.
  - It shows that under random allocation, a like-for-like comparator arm
    brings the test back to its correct error rate (Cell B).
  - It shows that targeted recruitment, especially extreme truncation, causes
    further problems, and compares how Oldham and the MLM hold up (Cells C
    and D).
  - A threshold sweep repeats Cells C and D at four recruitment cut-offs
    (top 25%, 5%, 1% and 0.2%). This shows how Type 1 error and power change
    as recruitment becomes more severe.
  - A confirmatory diagnostic compares the observed variance of the Fisher z
    contrast with the variance the theory predicts. This explains why the
    comparator test goes wrong under truncation.
  - An optional appendix repeats the analysis with a different definition of
    true_ide (one that shrinks only the variation not explained by baseline).
    This checks that the conclusions do not depend on how the IDE is defined.


OUTPUTS
--------------------------------------------------------------------------------

```
  Data/K_<K>/IDE results tables.xlsx, with these sheets:
      Table 2         Type 1 error and power, Oldham vs. MLM, Cells A to D
      Table A1        Power under the alternative true_ide definition
                      (appendix; left out if skip_appendix = TRUE)
      Table A2        Empirical vs. nominal Fisher z variance
      Figure A1       Threshold sweep results behind the severity plot
      MLM exclusions  Replications with no usable p-value, or where the
                      model did not converge

  Plots/          The severity plot: Type 1 error and power against
                  recruitment threshold, Oldham vs. MLM, Cells C and D.

  Data/K_<K>/     Cached .rds files (replication results). Their file
                  names carry a hash of the parameters used, so an old cache
                  is never reloaded by mistake after a parameter changes.
```


FILES IN THIS REPOSITORY
--------------------------------------------------------------------------------

```
  2026-09-28 IDE Simulations.R   The current simulation script
  RunScript.R                    Sources the simulation script
  README.md                      This file
```

The script follows a Pascal-like structure:

```
  PART 1  Setup and global parameters
  PART 2  Functions library (definitions only)
  PART 3  Functions verification
  PART 4  Simulation, cache and diagnostics
  PART 5  Confirmatory diagnostic (Fisher z reference variance)
  PART 6  Figures and tables
```


RUNNING THE SCRIPT
--------------------------------------------------------------------------------

Open the RStudio project, whose name must match the working directory, and
source RunScript.R. The script installs any missing packages itself: MASS,
ggplot2, parallel, compiler, psych, cocor, openxlsx, grid and lavaan.

The main settings are in PART 1, under RUN-CONTROL FLAGS. You can also set
them in your session before sourcing the script:

```
  K               Number of replications (default 1e6)
  sweep_K_max     Upper limit on replications for the threshold sweep
                  (default 1e5)
  skip_appendix   TRUE skips the costly appendix pass
  force_rerun_*   TRUE rebuilds that cache (OLDHAM, MLM, SWEEP)
  cache_version   Change this to invalidate every existing cache
```

A full run from scratch takes many hours, because the MLM fits are slow. For
a quick end-to-end test, set K to a small value (e.g. K <- 100) before
sourcing. The master copy of the script and its Data/ cache of results are
kept on the author's own machine. Data/ is not part of this repository.
