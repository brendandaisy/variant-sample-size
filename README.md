Sample size calculations for tracking genetic variants
================
Bren Case

## Calculator inputs

- `tmax` The number of time points (days, weeks, etc.) to do
  surveillance
- `q0` The initial proportion of the variant when sampling begins. For
  example, for a novel variant of concern a reasonable setting would be
  $1/N$ where $N$ is the current number of confirmed hospitalizations in
  your health care system
- `r` The true growth advantage. Because the growth advantage is the
  thing we are trying to estimate, we of course don’t know what it is
  (isn’t frequentist statistics fun?). Instead, there are two good ways
  to think about this input
  1.  It is the fastest growing variant you would like to learn with a
      given accuracy. For example, setting `r` to an Omicron-like growth
      (e.g. $0.84/week$) means that the sample size will have the
      intended accuracy for all growth advantages between $0$ and $0.84$
      in the future
  2.  A best “first guess” for the growth advantage. For example, you
      may first collect sequences to accurately detect an Omicron-like
      growth rate, however, after a few weeks it becomes clear the VOC
      is actually growing more slowly. Reducing the input for `r` and
      rerunning the calculator would then reduce the necessary sample
      size for the coming weeks.

## Examples

First we’ll need to load the following scripts:

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

Then, we set the following inputs. We would like to detect a VOC growth
advantage within `tmax=8` weeks, with initial proportion `q0=1/10,000`.
For the nominal growth advantage, we use an Omicron-like estimate of
$r=0.12$ per day, and convert it to a weekly timescale.

``` r
tmax <- 8
q0 <- 1e-4
r <- 0.12*7
```

Then, the optimal number of sequences to collect for an accuracy of
$15\%$ and with the default $95\%$ confidence, we call the function

``` r
nopt_single_param(tmax, q0, r, error=0.15*r)
```

    ## [1] 34.95864

If instead we only wanted $25\%$ and with $85\%$ confidence, use

``` r
nopt_single_param(tmax, q0, r, error=0.25*r, level=0.85)
```

    ## [1] 6.788958

As you can see, lowering our tolerance for accuracy or confidence level
can significantly reduce the number of sequences needed to collect each
week.

### Specifying a prevalence threshold

Instead of specifying *how long* we want to perform sampling, a more
practical way to think about things might be to say *how prevalent* we
are willing to let the VOC become before we need to know how fast it is
growing.

To do so, use the `find_tmax` function before calling the calculator.
For example, let’s say we want to estimate our growth advantage by the
time the VOC has reached 3% of circulating variants in the population:

``` r
qmax <- 0.03
q0 <- 1e-4
r <- 0.12*7

(tmax <- find_tmax(qmax, q0, r))
```

    ## [1] 7

``` r
nopt_single_param(tmax, q0, r, error=0.15*r)
```

    ## [1] 100.9203

As you can see, the VOC reaches 3% within 7 weeks, and 101 sequences are
needed per week in this case. If instead we are willing to wait until it
has reached 8%, we only need 14 weekly sequences:

``` r
tmax <- find_tmax(0.08, q0, r)
nopt_single_param(tmax, q0, r, error=0.15*r)
```

    ## [1] 13.54256
