# Extract and Reconstruct Wavelet Components from Significant Periods

Reconstructs time series components from the wavelet transform for
user-specified significant periods. Optionally includes a residual
component representing all non-significant scales, ensuring complete
signal decomposition.

## Usage

``` r
extract_wavelet_components(
  wave,
  signif_periods,
  scale,
  dj = 0.25,
  dt = 1,
  variable_sd,
  variable_mean = 0,
  Cdelta = 0.776,
  w0_0 = pi^(-1/4),
  include_residual = TRUE
)
```

## Arguments

- wave:

  Complex matrix. Wavelet transform coefficients with dimensions (scales
  x time). Output from continuous wavelet transform.

- signif_periods:

  Integer vector. Indices of significant period scales to reconstruct.
  Each index corresponds to a row in wave.

- scale:

  Numeric vector. Scale values corresponding to rows of wave.

- dj:

  Numeric scalar. Scale resolution parameter used in wavelet transform
  (default = 0.25). Smaller values give finer resolution.

- dt:

  Numeric scalar. Time step of the original series (default = 1).

- variable_sd:

  Numeric scalar. Standard deviation of the original time series before
  standardization. Used to restore original units.

- variable_mean:

  Numeric scalar. Mean of the original time series before
  standardization (default = 0). Added back to noise component to ensure
  complete reconstruction: rowSums(output) = original_signal.

- Cdelta:

  Numeric scalar. Reconstruction constant for Morlet wavelet (default =
  0.776). Depends on wavelet type; see Torrence & Compo (1998) Table 2.

- w0_0:

  Numeric scalar. Normalization constant psi_0(0) = pi^(-1/4) for Morlet
  wavelet (default = pi^(-1/4)).

- include_residual:

  Logical. If TRUE (default), includes a "Noise" component representing
  the sum of all non-significant scales, ensuring complete signal
  decomposition where rowSums(output) equals the reconstructed signal.

## Value

A matrix with columns: Component_1, Component_2, ...: One column per
significant period Noise: (if include_residual = TRUE) Sum of all
non-significant scales Number of rows equals the length of the original
time series. Property: rowSums(output) approximately equals the original
signal.

## Details

The reconstruction formula (Torrence & Compo 1998, Eq. 11) is:

x_n = (delta_j \* sqrt(delta_t)) / (C_delta \* psi_0(0)) \* sum over j
of Re(W_n(s_j)) / sqrt(s_j)

This function performs a complete wavelet decomposition by: 1.
Extracting each significant period as a separate component 2. Summing
all remaining (non-significant) scales into a "Noise" component 3.
Ensuring: original_signal = sum(significant_components) + Noise

This differs from partial reconstruction, which would only sum
significant scales and lose information from non-significant scales.

## References

Torrence, C. and Compo, G.P. (1998). A Practical Guide to Wavelet
Analysis. Bulletin of the American Meteorological Society, 79(1), 61-78.
See Section 3i and Table 2.

## Examples

``` r
if (FALSE) { # \dontrun{
# After performing wavelet analysis
components <- extract_wavelet_components(
  wave = wave_transform,
  signif_periods = c(5, 10, 20),  # indices of significant scales
  scale = scale_vector,
  dj = 0.25,
  dt = 1,
  variable_sd = sd(original_data),
  variable_mean = mean(original_data),
  Cdelta = 0.776,
  include_residual = TRUE
)

# Verify complete reconstruction
reconstructed <- rowSums(components)
plot(original_data)
lines(reconstructed, col = "red")

# Components include:
# - Component_1: First significant period
# - Component_2: Second significant period
# - Component_3: Third significant period
# - Noise: All non-significant scales combined + mean
} # }
```
