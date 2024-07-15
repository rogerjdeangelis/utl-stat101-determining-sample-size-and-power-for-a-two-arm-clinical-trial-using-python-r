# utl-stat101-determining-sample-size-and-power-for-a-two-arm-clinical-trial-using-python-r
Stat101 determining sample size and power for a two arm clinical trial using python r
    %let pgm=utl-stat101-determining-sample-size-and-power-for-a-two-arm-clinical-trial-using-python-r;

    Stat101 determining sample size and power for a two arm clinical trial using python r

    github
    https://tinyurl.com/mvs5m6aa
    https://github.com/rogerjdeangelis/utl-stat101-determining-sample-size-and-power-for-a-two-arm-clinical-trial-using-python-r

    Consider a two arm trial of aspirin vs placebo

    Consider a two arm clinical trial with equal sample size and common variance.
    The means of arm1 and arm2 may be different but the two standard deviations are the same.
    This is the simplest case

      Two Solution

         1 hand coded n on ttest
         2 r pwr package

    It is interesting that the t-test is stronly related to this problem and many
    powerful multivariate statistical procedures like manova, mancova and
    Hotelling's T-squared. The tetest handles dimensionality 2.

     Study Design Parameters

      1. Statistical Power:

          We aim for an 80% probability of detecting a statistically significant
          difference or effect when such a difference truly exists.

     2.  Expected Effect Size:

          - We anticipate a Cohen's d of 0.5.
          - If the standard deviation is 10 for both arms
            and the mean difference is 5 then the
            Cohen's d standardized difference is 5/10 or .5
            This eliminates the need for standard deviations

     3.  Significance Level:

          We want a two sided significance level of 5%  (alpha = 0.05).

    For more information on power, search this in perplexity
    How is statistical power related to the non-centrality of the t distribution

    Also note there are several packages in R for more complex statistical designs

    /*   _                     _                 _          _          _   _            _
    / | | |__   __ _ _ __   __| |   ___ ___   __| | ___  __| |  _ __  | |_| |_ ___  ___| |_
    | | | `_ \ / _` | `_ \ / _` |  / __/ _ \ / _` |/ _ \/ _` | | `_ \ | __| __/ _ \/ __| __|
    | | | | | | (_| | | | | (_| | | (_| (_) | (_| |  __/ (_| | | | | || |_| ||  __/\__ \ |_
    |_| |_| |_|\__,_|_| |_|\__,_|  \___\___/ \__,_|\___|\__,_| |_| |_| \__|\__\___||___/\__|

    */


    Consider a two arm clinical trial with equal sample size and common variance.
    The means of arm1 and arm2 may be different but the two standard deviations are the same.


                      s         +  s
                       placebo      aspirin
      s          =   ----------------------
       poooled                 2


      Definition of ttset (stat101)

                        mean      -  mean
      t                  placebo      aspirin
       statistic    =   ----------------------
                          .___________________
                          |
                          |  2       /    \
                      |\  | s       |  2  |
                        \ |  pooled |  -  |
                         \|          \ n  /



                  mean      -  mean
                   placebo      aspirin
      Let d =     ----------------------
                         .________
                         |
                         |  2
                     |\  | s
                       \ |  pooled
                        \|
      Thus

      t                      d
       statistic   =   --------------
                           .________
                           |
                           |  /    \
                       |\  | |  2  |
                         \ | |  -  |
                          \|  \ n  /


       Let set the significance of alpha and the power of beta


       t_alpha  = stats.t.ppf(1 - alpha/2, n)
       t_beta   = stats.t.ppf(power, n)       (80% chanc we found a sinificant difference when thee is one)

       Unlike the normal distribution the t distribution critical vales depend
       on the sample size, so we will need to iterate.


       Basicall we need to solve for n by iterating

                                   d
                             --------------
       t_alpha + t_beta  =      .________
                                 |
                                 |  /    \
                             |\  | |  2  |
                               \ | |  -  |
                                \|  \ n  /


       see above for info in poer                        d
                                                    --------------
                                                       .________
       probt(1 - alpha/2, n) + probt((power, n) =       |
                                                        |  /    \
                                                    |\  | |  2  |
                                                      \ | |  -  |
                                                       \|  \ n  /


                                                        2
                                                 2   d * n
       (probt(1 - alpha/2, n) + probt((power, n))  = -----
                                                       2

       Problw is that n is on both sides.

       There may be a closed form solution for n using sympy but
       I suspect it would be very messy, so lets use an
       iterative solution.

       For large sample size the solution does not have n on
       left side(normal approx) of the equation and we can fomulate(with erf)
       a closed form solution or directly solve with sas probit(alpha)
       instead of probt(1-alpha/2,n)

    /* _                 _   _                           _
    (_) |_ ___ _ __ __ _| |_(_)_   _____    ___ ___   __| | ___
    | | __/ _ \ `__/ _` | __| \ \ / / _ \  / __/ _ \ / _` |/ _ \
    | | ||  __/ | | (_| | |_| |\ V /  __/ | (_| (_) | (_| |  __/
    |_|\__\___|_|  \__,_|\__|_| \_/ \___|  \___\___/ \__,_|\___|

    */
    %utl_pybegin;
    parmcards4;
    from scipy import stats
    import numpy as np
    def calculate_sample_size(alpha, power, effect_size):
        df = float('inf')  # Assume large sample size initially
        t_alpha = stats.t.ppf(1 - alpha/2, df)
        t_beta = stats.t.ppf(power, df)

        n = 2 * ((t_alpha + t_beta) / effect_size)**2

        # Round up to the nearest integer
        n = int(np.ceil(n))

        # Iterate to refine the estimate using the correct degrees of freedom
        for _ in range(5):  # Usually converges within a few iterations
            df = 2*n - 2
            t_alpha = stats.t.ppf(1 - alpha/2, df)
            t_beta = stats.t.ppf(power, df)
            n_new = int(np.ceil(2 * ((t_alpha + t_beta) / effect_size)**2))
            if n_new == n:
                break
            n = n_new

        return n

    # Example usage
    alpha = 0.05
    power = 0.8
    effect_size = 0.5

    sample_size = calculate_sample_size(alpha, power, effect_size)
    print(f"Required sample size per group: {sample_size}")
    ;;;;
    %utl_pyend;

    /**************************************************************************************************************************/
    /*                                                                                                                        */
    /*        alpha = 0.05                                                                                                    */
    /*        power = 0.8                                                                                                     */
    /*        effect_size = 0.5                                                                                               */
    /*                                                                                                                        */
    /*        Required sample size per group: 64                                                                              */
    /*                                                                                                                        */
    /**************************************************************************************************************************/

    /*___                                                 _
    |___ \   _ __   _ ____      ___ __   _ __   __ _  ___| | ____ _  __ _  ___
      __) | | `__| | `_ \ \ /\ / / `__| | `_ \ / _` |/ __| |/ / _` |/ _` |/ _ \
     / __/  | |    | |_) \ V  V /| |    | |_) | (_| | (__|   < (_| | (_| |  __/
    |_____| |_|    | .__/ \_/\_/ |_|    | .__/ \__,_|\___|_|\_\__,_|\__, |\___|
                   |_|                  |_|                         |___/
    */

    %utl_rbeginx;
    parmcards4;
    library(pwr)

    # Set parameters
    alpha <- 0.05  # Significance level
    power <- 0.80  # Desired power
    d <- 0.5       # Expected effect size (Cohens d)

    # Calculate sample size
    sample_size <- pwr.t.test(d = d,
                              sig.level = alpha,
                              power = power,
                              type = "two.sample",
                              alternative = "two.sided")

    str(sample_size);
    sample_size

    # Print required sample size per group
    print(ceiling(sample_size$n))

    # Check power for a given sample size
    n <- 64  # Sample size per group
    achieved_power <- pwr.t.test(n = n,
                                 d = d,
                                 sig.level = alpha,
                                 type = "two.sample",
                                 alternative = "two.sided")

    str(achieved_power)
    achieved_power
    # Print achieved power
    print(achieved_power$power)
    ;;;;
    %utl_rendx;


    /**************************************************************************************************************************/
    /*                                                                                                                        */
    /*      Two-sample t test power calculation                                                                               */
    /*                                                                                                                        */
    /*               n = 63.76561      ---SAMPLE SIZE PER ARM 128 TOTAL SAME AS ITERATIVE CODE ABOVE ---                      */
    /*                                                                                                                        */
    /*               d = 0.5                                                                                                  */
    /*       sig.level = 0.05                                                                                                 */
    /*           power = 0.8                                                                                                  */
    /*     alternative = two.sided                                                                                            */
    /*                                                                                                                        */
    /*    Using the n lets check the power                                                                                    */
    /*                                                                                                                        */
    /*         Two-sample t test power calculation                                                                            */
    /*                                                                                                                        */
    /*                  n = 64                                                                                                */
    /*                  d = 0.5                                                                                               */
    /*          sig.level = 0.05                                                                                              */
    /*                                                                                                                        */
    /*              power = 0.8014596  ---POWER---                                                                            */
    /*                                                                                                                        */
    /*        alternative = two.sided                                                                                         */
    /*                                                                                                                        */
    /**************************************************************************************************************************/

    /*              _
      ___ _ __   __| |
     / _ \ `_ \ / _` |
    |  __/ | | | (_| |
     \___|_| |_|\__,_|

    */

