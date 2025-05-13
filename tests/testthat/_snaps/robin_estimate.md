# robin_wt snapshot

    Code
      robin_wt(data = data_sim, estimand = list(tx_colname = tx_colname, comparison = c(
        "trt.1", "trt.2")), design = list(randomization_var_colnames = randomization_var_colnames,
        randomization_table = randomization_table), estimated_propensity = FALSE,
      outcome_model = list(formula = y ~ xb + xc, family = gaussian()))
    Output
      Method:  Inverse Probability Weighting 
      Model :  y ~ xb + xc 
      Family:  gaussian 
      Randomization Probabilities (among the entire concurrent and eligible (ECE) population): 
        t subtype trt.1 trt.2 Sample.Size Proportion
      1 1       0   0.5  0.50          47       0.09
      2 1       1   0.5  0.20         105       0.21
      3 2       0   0.5  0.50          46       0.09
      4 2       1   0.5  0.15         133       0.27
      5 3       0   0.5  0.50          11       0.02
      6 3       1   0.5  0.20         158       0.32
      
      Nominal Level:  0.05 
      ---------------------------
      Marginal Mean: 
            Estimate Std.Err Z Value lower.CL upper.CL
      trt.1    2.243   0.124    18.1    2.000    2.486
      trt.2    5.101   0.302    16.9    4.509    5.693
      
      ---------------------------
      Treatment Effect: 
      Contrast:  difference 
                    Estimate Std.Err Z Value lower.CL upper.CL Pr(>|z|)    
      trt.2 - trt.1    2.858   0.335    8.53    2.201    3.515   <2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

---

    Code
      robin_wt(data = data_sim, estimand = list(tx_colname = tx_colname, comparison = c(
        "trt.1", "trt.2")), design = list(randomization_var_colnames = randomization_var_colnames,
        randomization_table = NULL), estimated_propensity = TRUE, outcome_model = list(
        formula = y ~ xb + xc, family = gaussian()))
    Output
      Method:  Inverse Probability Weighting 
      Model :  y ~ xb + xc 
      Family:  gaussian 
      Estimated Propensity Score is used.
      Estimated Randomization Probabilities (among the entire concurrent and eligible (ECE) population): 
        t subtype trt.1 trt.2 Sample.Size Proportion
      1 1       0  0.53  0.47          47       0.09
      2 1       1  0.45  0.25         105       0.21
      3 2       0  0.50  0.50          46       0.09
      4 2       1  0.49  0.14         133       0.27
      5 3       0  0.64  0.36          11       0.02
      6 3       1  0.45  0.23         158       0.32
      
      Nominal Level:  0.05 
      ---------------------------
      Marginal Mean: 
            Estimate Std.Err Z Value lower.CL upper.CL
      trt.1    2.272   0.127    17.9    2.024    2.521
      trt.2    5.140   0.289    17.8    4.573    5.706
      
      ---------------------------
      Treatment Effect: 
      Contrast:  difference 
                    Estimate Std.Err Z Value lower.CL upper.CL Pr(>|z|)    
      trt.2 - trt.1    2.867   0.323    8.87    2.234    3.501   <2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# robin_ps snapshot

    Code
      robin_ps(data = data_sim, estimand = list(tx_colname = tx_colname, comparison = c(
        "trt.1", "trt.2")), design = list(randomization_var_colnames = randomization_var_colnames,
        randomization_table = NULL), stratify_by = "s12", outcome_model = list(
        formula = y ~ xb + xc, family = gaussian()))
    Condition
      Warning in `consistency_check()`:
      Assignment probabilities for the arms being compared must be constant within each level of stratify_by.
      This check is skipped because randomization_table was not provided.
    Output
      Method:  Post Stratification 
      Post stratification is done by variable s12 specified by stratify_by.
      Model :  y ~ xb + xc 
      Family:  gaussian 
      Estimated Randomization Probabilities (among the entire concurrent and eligible (ECE) population): 
        trt.1 trt.2 Sample.Size Proportion Stratum
      1  0.53  0.47         104       0.21       a
      2  0.45  0.24         263       0.53       b
      3  0.49  0.14         133       0.27       c
      
      Nominal Level:  0.05 
      ---------------------------
      Marginal Mean: 
            Estimate Std.Err Z Value lower.CL upper.CL
      trt.1    2.273   0.125    18.2    2.029    2.517
      trt.2    5.131   0.292    17.5    4.557    5.704
      
      ---------------------------
      Treatment Effect: 
      Contrast:  difference 
                    Estimate Std.Err Z Value lower.CL upper.CL Pr(>|z|)    
      trt.2 - trt.1    2.858   0.324    8.82    2.223    3.493   <2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# robin_ps snapshot2

    Code
      robin_ps(data = data_sim, estimand = list(tx_colname = tx_colname, comparison = c(
        "trt.1", "trt.2")), design = list(randomization_var_colnames = randomization_var_colnames,
        randomization_table = randomization_table), stratify_by = NULL,
      outcome_model = list(formula = y ~ xb + xc, family = gaussian()))
    Output
      Method:  Post Stratification 
      Post stratification is done by the joint levels of the randomization variables specified by randomization_var_colnames.
      Model :  y ~ xb + xc 
      Family:  gaussian 
      Randomization Probabilities (among the entire concurrent and eligible (ECE) population): 
        t subtype trt.1 trt.2 Sample.Size Proportion
      1 1       0   0.5  0.50          47       0.09
      2 1       1   0.5  0.20         105       0.21
      3 2       0   0.5  0.50          46       0.09
      4 2       1   0.5  0.15         133       0.27
      5 3       0   0.5  0.50          11       0.02
      6 3       1   0.5  0.20         158       0.32
      
      Nominal Level:  0.05 
      ---------------------------
      Marginal Mean: 
            Estimate Std.Err Z Value lower.CL upper.CL
      trt.1    2.273   0.125    18.2    2.029    2.517
      trt.2    5.131   0.292    17.5    4.557    5.704
      
      ---------------------------
      Treatment Effect: 
      Contrast:  difference 
                    Estimate Std.Err Z Value lower.CL upper.CL Pr(>|z|)    
      trt.2 - trt.1    2.858   0.324    8.82    2.223    3.493   <2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

