# Test Whether an Object Is a Finite Integer Scalar

Checks whether an object is a numeric scalar representing a finite
integer value. This helper is intended for lightweight input validation
where strict integer typing is not required but integer-valued numerics
are acceptable.

## Usage

``` r
.is_int_scalar(x)

.is_int_scalar(x)
```

## Arguments

- x:

  Object.

## Value

Logical scalar indicating whether `x` is a finite integer-valued scalar.

Logical scalar.

## Details

The function returns `TRUE` if and only if:

- `x` is numeric,

- `x` has length 1,

- `x` is finite (not `NA`, `NaN`, or `Inf`),

- `x` has no fractional component (`x %% 1 == 0`).

Logical, character, and integer vectors of length greater than one will
return `FALSE`.
