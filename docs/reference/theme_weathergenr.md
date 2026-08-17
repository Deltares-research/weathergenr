# weathergenr house theme

A [`theme_light`](https://ggplot2.tidyverse.org/reference/ggtheme.html)
wrapper applied to every figure this package draws, so the PNGs a run
writes share one look.

[`evaluate_weather_generator()`](https://deltares-research.github.io/weathergenr/reference/evaluate_weather_generator.md)
returns its diagnostic plots as ggplot objects, so this is exported:
without it, a caller re-rendering one of those plots cannot reproduce
the package's styling.

## Usage

``` r
theme_weathergenr(base_size = 12, base_family = "")
```

## Arguments

- base_size:

  Numeric. Base font size in points. Default 12.

- base_family:

  Character. Base font family. Default `""`, the device default.

## Value

A ggplot2 theme object.

## Details

`base_size` is 12 rather than
[`theme_light()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)'s
native 11 because 12 is what the evaluation pipeline has always
produced. The title and subtitle sizes are
[`rel`](https://ggplot2.tidyverse.org/reference/element.html) rather
than absolute points so they track `base_size`; they were previously
fixed at 14 and 10 and would have been stranded by any change to it.

Legend text is deliberately not set here. Exactly one exported figure
has a visible legend, so its styling stays local to that plot rather
than becoming a package-wide rule inferred from a single case.

## Examples

``` r
library(ggplot2)
#> Warning: package 'ggplot2' was built under R version 4.6.1
ggplot(mtcars, aes(wt, mpg)) +
  geom_point() +
  theme_weathergenr()

```
