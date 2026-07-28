%let pgm=utl-stat101-determining-sample-size-and-power-for-a-two-arm-clinical-trial-using-python-r;

/* -------------------------------------------------------------------------
   Two-arm clinical trial: sample size and power (aspirin vs placebo).

   Native-SAS port of the calculate_sample_size routine in
   utl-stat101-...-using-python-r.sas.  The original spells the algorithm
   out in a %utl_pybegin Python block (scipy.stats.t.ppf) and re-checks it
   with the R pwr package.  Here the same iterative t-distribution solve is
   expressed with SAS's own TINV / PROBIT functions -- the direct analogs of
   scipy.stats.t.ppf -- keeping the repo's variable names (alpha, power,
   effect_size, t_alpha, t_beta, n) and its 5-iteration refinement loop.

   Study design (from the repo):
     alpha       = 0.05   two-sided significance
     power       = 0.80
     effect_size = 0.5    Cohen's d
   Expected result: n = 64 per arm (128 total), matching the repo's
   documented Python output and the R pwr.t.test n = 63.76561 -> ceil 64.
   ------------------------------------------------------------------------- */

data sample_size;

  alpha       = 0.05;   /* two-sided significance level          */
  power       = 0.80;   /* desired power                         */
  effect_size = 0.5;    /* Cohen's d = mean diff / pooled sd     */

  /* Initial estimate: normal approximation (df -> infinity),
     i.e. probit() stands in for t.ppf at large sample size.     */
  t_alpha = probit(1 - alpha/2);
  t_beta  = probit(power);
  n = ceil( 2 * ((t_alpha + t_beta) / effect_size)**2 );

  /* Refine with the exact t quantiles, df = 2n - 2, iterating as
     in the repo's Python loop (t critical values depend on n).  */
  do _iter = 1 to 5;
    df      = 2*n - 2;
    t_alpha = tinv(1 - alpha/2, df);
    t_beta  = tinv(power, df);
    n_new   = ceil( 2 * ((t_alpha + t_beta) / effect_size)**2 );
    if n_new = n then leave;
    n = n_new;
  end;

  put "Required sample size per group: " n;
run;

proc print data=sample_size noobs;
  var alpha power effect_size n;
run;
