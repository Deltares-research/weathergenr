# print and summary methods are stable

    Code
      print(out)
    Output
      
      === Precipitation QM Diagnostics ===
      
      Mean and variance changes:
       Statistic Original Adjusted Target Achieved Error
            Mean    10.36    12.43   +20%     +20%   +0%
        Variance    49.24    70.90   +44%     +44%   +0%
      
      Wet/dry day counts:
       Category Original Adjusted Difference
       Dry days       36       36          0
       Wet days       84       84          0
      
      Percentiles (wet days):
         Percentile Original Adjusted Change
       P50 (Median)     9.00    10.80   +20%
                P90    20.00    24.00   +20%
                P95    24.00    28.80   +20%
                P99    28.68    34.42   +20%
      
      Extreme tails:
                   Tail Threshold change Mean change
          Top 10% (P90)             +20%        +20%
           Top 5% (P95)             +20%        +20%
           Top 1% (P99)             +20%        +20%
       Top 0.1% (P99.9)             +20%        +20%
      
      Spell lengths:
             Type Original (days) Adjusted (days) Change
       Dry spells             1.0             1.0    +0%
       Wet spells             2.3             2.3    +0%
      
      Monthly mean changes:
       Month Target Achieved Error
         Jan   +20%     +20%   +0%
         Feb   +20%     +20%   +0%
         Mar   +20%     +20%   +0%
         Apr   +20%     +20%   +0%
         May   +20%     +20%   +0%
         Jun   +20%     +20%   +0%
         Jul   +20%     +20%   +0%
         Aug   +20%     +20%   +0%
         Sep   +20%     +20%   +0%
         Oct   +20%     +20%   +0%
         Nov   +20%     +20%   +0%
         Dec   +20%     +20%   +0%
      

---

    Code
      summary(out)
    Output
      
      === Precipitation QM Summary ===
      
      Targets vs achieved:
        Mean:     achieved 1.200 | target 1.200 | error +0.000
        Variance: achieved 1.440 | target 1.440 | error +0.000
      
      Wet/dry day preservation:
        Dry days change: 0
        Wet days change: 0
      
      Quantile changes:
        Max absolute change: 20.0%
      
      Extreme tail ratios (mean across thresholds):
        Threshold: 1.200
        Tail mean: 1.200
      

