# print and summary methods are stable

    Code
      print(out)
    Output
      
      === Precipitation QM Diagnostics ===
      
      Mean and variance changes:
       Statistic Original Adjusted Target Achieved Error
            Mean    10.36    12.95   +25%     +25%   +0%
        Variance    49.24    76.93   +56%     +56%   +0%
      
      Wet/dry day counts:
       Category Original Adjusted Difference
       Dry days       36       36          0
       Wet days       84       84          0
      
      Percentiles (wet days):
         Percentile Original Adjusted Change
       P50 (Median)     9.00    11.25   +25%
                P90    20.00    25.00   +25%
                P95    24.00    30.00   +25%
                P99    28.68    35.85   +25%
      
      Extreme tails:
                   Tail Threshold change Mean change
          Top 10% (P90)             +25%        +25%
           Top 5% (P95)             +25%        +25%
           Top 1% (P99)             +25%        +25%
       Top 0.1% (P99.9)             +25%        +25%
      
      Spell lengths:
             Type Original (days) Adjusted (days) Change
       Dry spells             1.0             1.0    +0%
       Wet spells             2.3             2.3    +0%
      
      Monthly mean changes:
       Month Target Achieved Error
         Jan   +25%     +25%   +0%
         Feb   +25%     +25%   +0%
         Mar   +25%     +25%   +0%
         Apr   +25%     +25%   +0%
         May   +25%     +25%   +0%
         Jun   +25%     +25%   +0%
         Jul   +25%     +25%   +0%
         Aug   +25%     +25%   +0%
         Sep   +25%     +25%   +0%
         Oct   +25%     +25%   +0%
         Nov   +25%     +25%   +0%
         Dec   +25%     +25%   +0%
      

---

    Code
      summary(out)
    Output
      
      === Precipitation QM Summary ===
      
      Targets vs achieved:
        Mean:     achieved 1.250 | target 1.250 | error +0.000
        Variance: achieved 1.563 | target 1.562 | error +0.000
      
      Wet/dry day preservation:
        Dry days change: 0
        Wet days change: 0
      
      Quantile changes:
        Max absolute change: 25.0%
      
      Extreme tail ratios (mean across thresholds):
        Threshold: 1.250
        Tail mean: 1.250
      

