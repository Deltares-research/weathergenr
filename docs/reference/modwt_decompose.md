# Maximal Overlap Discrete Wavelet Transform (MODWT)

Performs the Maximal Overlap Discrete Wavelet Transform using a
specified wavelet filter. MODWT is shift-invariant and provides additive
multi-resolution decomposition at dyadic scales (powers of 2).

This function wraps
[`waveslim::modwt()`](https://rdrr.io/pkg/waveslim/man/modwt.html) with
input validation and a standardized return object.

## Usage

``` r
modwt_decompose(
  x,
  filter = c("la8", "haar", "d4", "d6", "d8", "la16"),
  n_levels = NULL,
  boundary = "periodic"
)
```

## Arguments

- x:

  Numeric vector. Input time series. Must be regularly spaced and
  contain no missing values.

- filter:

  Character. Wavelet filter name. One of `"la8"` (default), `"haar"`,
  `"d4"`, `"d6"`, `"d8"`, `"la16"`.

- n_levels:

  Integer scalar or NULL. Number of decomposition levels (J). If NULL, a
  conservative default is used based on series length and filter width.

- boundary:

  Character. Boundary handling method. Only `"periodic"` is supported.

## Value

A list with class `"modwt_result"`. The object contains: -
`coefficients`: coefficients returned by
[`waveslim::modwt()`](https://rdrr.io/pkg/waveslim/man/modwt.html)
(W1...WJ and VJ) - `filter`: filter name used - `n_levels`: number of
decomposition levels (J) - `n`: series length - `boundary`: boundary
method - `filter_length`: length of the wavelet filter

## Details

The recommended maximum number of levels is capped using a stability
heuristic based on series length and filter width. The cap is intended
to reduce boundary artifacts and unstable behavior on short records; it
is not a mathematical requirement of MODWT.
