# robin_wt snapshot

    Code
      robin_wt(y ~ xb + xc, data = example, treatment = "treatment", prob_mat = prob_mat,
      treatments_for_compare = c("1", "3"), contrast = "difference")
    Output
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

