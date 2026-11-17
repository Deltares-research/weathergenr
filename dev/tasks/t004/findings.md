# Findings: interannual variability under `enforce_target_mean`

Task `t004` (formerly `t2611161814`), 2026-11-16. Source of every number: the full run of
`dev/scripts/perturbation-figures.R` (20 realizations x 30 calendar years, 2030-2059, Ntoum
fixture, six cells; generator seed 100, perturbation seed 42). Raw tables are in
`output/perturbation-figures/` (git-ignored): `phase0-variance.csv`, `phase0-summary.csv`,
`wet-day-fractions.csv`, `checks.md`, `run-info.txt`.

Validation status: real-data exercised (fixture), not scientifically qualified.

## Conclusion

1. The suspected defect is real and material. With `enforce_target_mean = TRUE` (the default) the
   interannual CV of annual precipitation falls to about 1% of its baseline value. Median CV ratio:
   0.012, range 0.001-0.028. The falsifier threshold was 0.9.
2. The cause is the per-(year, month) rescale, not the mapping. With enforcement off, the mapping is an
   exact 1.15 x scaling, so the CV ratio is 1.000 and the anomaly correlation is 1.000.
3. The same rescale also compresses daily extremes and shifts daily structure, so the defect reaches
   past annual totals (see below). Whether and how to change it is an owner decision.

## Phase 0 table

Scenario: S2 hot-wet in step mode (+4 °C, wet-day mean x1.15, variance scaled with the mean). The
table summarizes 120 realization x cell pairs; ratios are perturbed over baseline within the same pair.
Wet-day mean uses the mapping's wet-day definition, > 0 mm.

| Metric | Baseline | Enforce TRUE: median (range) | Enforce FALSE: median (range) |
|---|---|---|---|
| Interannual CV, annual totals | 0.102 (0.066-0.148) | 0.0013 (0.0001-0.0035) | 0.102 (0.066-0.148) |
| **CV ratio** | 1 | **0.012 (0.001-0.028)** | **1.000 (1.000-1.000)** |
| SD ratio, annual totals | 1 | 0.014 (0.001-0.032) | 1.150 |
| SD across years of monthly wet-day mean (mm/day), mean over months | 1.68 (1.32-2.29) | 7e-16 (zero to rounding) | 1.94 (1.52-2.63) |
| Lag-1 autocorrelation, annual totals | 0.130 (-0.29 to 0.38) | -0.019 (-0.28 to 0.34) | 0.130 (identical) |
| Correlation of baseline and perturbed annual anomalies | 1 | -0.014 (-0.56 to 0.57) | 1.000 |

Under enforcement, the perturbed annual series carries almost no variance and no relation to the
baseline sequence of wet and dry years. Its lag-1 value and anomaly correlation are noise about zero.

## Mechanism

With `scale_var_with_mean = TRUE`, the target Gamma keeps the base shape and multiplies the scale by
f, so `qgamma(pgamma(x))` equals f x exactly. The enforcement step then multiplies each (year, month)
group by `clim_mean(m) * f / mean_mapped(y, m)`. This equals `clim_mean(m) / mean(y, m)`, the inverse
of that year's anomaly. Wet months are pulled down and dry months up, so:

- Each month-year's wet-day mean becomes exactly the climatological mean x f; its SD across years is
  zero to rounding.
- Monthly and annual totals then vary only through wet-day counts. On this fixture the mapping's
  threshold of 0 mm makes 99.8% of days wet, so that channel is essentially closed.
- The period total in each month is preserved exactly, because the climatological mean is the pooled
  mean. That is why mean-based checks pass (figure 6(a) tracks the factor exactly) while variability
  vanishes.

The fixture is close to the worst case: a climate with more dry days would keep more interannual
variability through wet-day counts. I did not run extra scenarios to measure that.

## Wet-day fractions

| Definition | Fixture observations, mean (per cell) | Baseline ensemble, mean |
|---|---|---|
| > 0 mm (mapping) | 0.998 (0.996-0.999) | 0.998 |
| >= 1 mm (statistics) | 0.815 (0.790-0.870) | 0.811 |

## Side effects of the rescale on daily structure (S2 step against S0)

Medians over cell x realization (structure statistics use >= 1 mm wet days).

| Statistic | S0 | Enforce TRUE | Enforce FALSE |
|---|---|---|---|
| Lag-1 autocorrelation, daily anomalies | 0.173 | 0.105 | 0.173 |
| Cross-cell daily correlation | 0.904 | 0.917 | 0.904 |
| P-T anomaly correlation | -0.308 | -0.307 | -0.308 |
| Dry-spell length (days) | 3.19 | 2.44 | 3.01 |
| Wet-spell length (days) | 12.98 | 13.31 | 13.66 |
| Interannual CV, annual P | 0.102 | 0.0013 | 0.102 |
| Rx1day change from S0 | — | -2.1% | +15.0% |

- **Persistence:** lag-1 drops because removing month-year mean anomalies removes persistent
  low-frequency variance.
- **Spells:** the spell and wet-fraction shifts with enforcement off are the expected effect of x1.15
  moving days across the 1 mm statistics threshold. The extra dry-spell shortening with it on comes
  from dry months being scaled up.
- **Figure 6, (ΔT = 0, factor 1.0) node:** this should be a no-op. Annual precipitation, PET and
  aridity change by 0.0000%, but **Rx1day changes by -14.9%** (per realization -18.5 to -11.1%).
- **Figure 3 (cell 1):** the M1 paired ratio is not flat near 1.15. It reads 1.10-1.15 for
  exceedance 0.5 to 0.05, 1.02 at 0.01 and 0.85 at 0.001. The same precipitation with enforcement
  off reads 1.128-1.145 for exceedance 0.5 to 0.05, 1.143 at 0.01 and 1.150 at 0.001. The small dip
  at the low end is pooling at >= 1 mm; the tail collapse is enforcement.
- **Figure 1(d):** realized annual ratios scatter about the specified path because perturbed totals
  are nearly constant while baseline totals vary; the perturbation does not add noise.

## Sanity checks (full run)

| Check | Result | Values |
|---|---|---|
| Fig 2: perturbed minus baseline monthly temperature | pass | S1 2.000000 °C, S2 4.000000 °C at every month and realization |
| Fig 3: M1 ratio flat near 1.15 below the tail | **fail** | 1.099-1.152 for exceedance >= 0.05, then 1.024 at 0.01 and 0.852 at 0.001 |
| Fig 3: M3 departs upward from 1.15 above exceedance 0.05 | **fail** | 1.144 at 0.02 and 0.01, 0.975 at 0.001; it stays above M1 (0.85) in the tail but never above 1.15 |
| Fig 6: (0, 1.0) node 0% in every panel | **fail (Rx1day)** | precip 0, PET 0, aridity 0, Rx1day -14.86% |
| Fig 6: annual precipitation tracks the factor | pass | -30.000% to +30.000% at factors 0.7-1.3, exact |
| Separable-grid equivalence | pass | full call against split runs: max difference 0 mm (precip and PET) |
| Self-check: CV ratio with enforcement off | pass | 1.000, sd ratio 1.150 |

Both fails are consequences of the enforcement step; neither was tuned.

## Other observations for the owner

- **M3 tail transform is discontinuous at the threshold.** `u' = 1 - (1 - u)^k` with `k = 1.2` maps
  0.95 to 0.9725, so mapped values jump at `extreme_prob_threshold` (the step in figure 1(b)). It is
  continuous only for `k = 1`. Combined with enforcement, M3's net tail effect on cell 1 is a ratio
  below 1 at exceedance 0.001.
- **Two wet-day definitions coexist.** The mapping uses > 0 mm (99.8% of days wet), while the
  statistics use >= 1 mm. The figure 3 line at exceedance 0.05 marks the tail threshold on the fitted
  Gamma of > 0 mm days, so it sits only approximately at 0.05 on the >= 1 mm axis.
- **The Gamma fits the > 0 mm sample poorly in some months.** Figure 1(a) shows November, cell 1: the
  empirical CDF departs from the method-of-moments fit, driven by the drizzle days.

## Owner decisions

1. **The `enforce_target_mean` behaviour.**
   - Options: rescale against each month-year's own mapped input mean times the factor, or enforce
     only over the pooled period.
   - Alternatively, default it to FALSE when `exaggerate_extremes = FALSE`, since the mapping already
     hits the mean exactly in that case.
   - Any change moves numeric output and needs the baseline gate.
2. **Whether the tail transform should be made continuous at `extreme_prob_threshold`.**
3. **Whether the default mapping threshold of 0 mm is appropriate** for daily data with near-daily
   drizzle.
4. **Acceptance of the figures**, and whether figure 4b belongs in the article set.
