# robin_ps snapshot

    Code
      robin_ps(y ~ xb + xc, data = example, treatment = "treatment", prob_mat = prob_mat,
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

