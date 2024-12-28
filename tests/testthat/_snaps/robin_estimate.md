# robin_wt snapshot

    Code
      robin_wt(y ~ xb + xc, data = example, treatment = "treatment", prob_mat = prob_mat,
      treatments_for_compare = c("1", "3"), contrast = "difference")
    Output
      Method:  Inverse Probability Weighting 
      Model :  y ~ xb + xc 
      Family:  gaussian 
      Randomization Probabilities (among the entire concurrent and eligible (ECE) samples): 
        Total Sample Size:  238 
        Unique.Level Sample.Size Proportion
      1 "0.5, 0.15"  "133"       "0.56"    
      2 "0.5, 0.3"   "105"       "0.44"    
      
      Nominal Level:  0.05 
      ---------------------------
      Marginal Mean: 
        Estimate Std.Err Z Value lower.CL upper.CL Pr(>|z|)    
      1    3.029   0.163  18.550    2.709     3.35   <2e-16 ***
      3    4.147   0.297  13.980    3.566     4.73   <2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      ---------------------------
      Treatment Effect: 
      Contrast:  difference 
            Estimate Std.Err Z Value lower.CL upper.CL Pr(>|z|)    
      3 - 1    1.119   0.326   3.430    0.480     1.76    6e-04 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# robin_ps snapshot

    Code
      robin_ps(y ~ xb + xc, data = example, treatment = "treatment", prob_mat = prob_mat,
      stratification = NULL, treatments_for_compare = c("1", "3"), contrast = "difference")
    Output
      Method:  Post Stratification 
      Model :  y ~ xb + xc 
      Family:  gaussian 
      Randomization Probabilities (among the entire concurrent and eligible (ECE) samples): 
        Total Sample Size:  238 
        Unique.Level Sample.Size Proportion
      1 "0.5, 0.15"  "133"       "0.56"    
      2 "0.5, 0.3"   "105"       "0.44"    
      
      Nominal Level:  0.05 
      ---------------------------
      Marginal Mean: 
        Estimate Std.Err Z Value lower.CL upper.CL Pr(>|z|)    
      1    3.023   0.171  17.670    2.687     3.36   <2e-16 ***
      3    4.149   0.270  15.340    3.619     4.68   <2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      ---------------------------
      Treatment Effect: 
      Contrast:  difference 
            Estimate Std.Err Z Value lower.CL upper.CL Pr(>|z|)    
      3 - 1    1.126   0.307   3.660    0.524     1.73  0.00025 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

---

    Code
      robin_ps(y ~ xb + xc, data = example, treatment = "treatment", prob_mat = NULL,
      stratification = "s12.2", treatments_for_compare = c("1", "2"), contrast = "difference")
    Condition
      Warning in `prob_strata_check()`:
      Prob is not provided. The result is not guaranteed to be valid.
    Output
      Method:  Post Stratification 
      Model :  y ~ xb + xc 
      Family:  gaussian 
      Randomization Probabilities (among the entire concurrent and eligible (ECE) samples): 
        Total Sample Size:  500 
        Unique.Level   Sample.Size Proportion
      1 "0.143, 0.489" "133"       "0.27"    
      2 "0.169, 0.446" "213"       "0.43"    
      3 "0.471, 0.529" "104"       "0.21"    
      4 "0.52, 0.46"   "50"        "0.1"     
      
      Nominal Level:  0.05 
      ---------------------------
      Marginal Mean: 
        Estimate Std.Err Z Value lower.CL upper.CL Pr(>|z|)    
      1    2.273   0.127  17.900    2.024     2.52   <2e-16 ***
      2    5.242   0.293  17.870    4.667     5.82   <2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      ---------------------------
      Treatment Effect: 
      Contrast:  difference 
            Estimate Std.Err Z Value lower.CL upper.CL Pr(>|z|)    
      2 - 1    2.969   0.356   8.330    2.270     3.67   <2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

