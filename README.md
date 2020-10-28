`gratis` <img src="man/figures/logo.PNG" align="right" height="210"/>
========

<!-- badges: start -->
  [![R build status](https://github.com/ykang/gratis/workflows/R-CMD-check/badge.svg)](https://github.com/ykang/gratis/actions)

  [![](https://cranlogs.r-pkg.org/badges/gratis)](https://CRAN.R-project.org/package=gratis)

  <!-- badges: end -->

The R package `gratis` (previously known as `tsgeneration`) provides efficient algorithms for generating time series with
diverse and controllable characteristics.

Installation
------------


### [CRAN version]( https://CRAN.R-project.org/package=gratis)


```r
install.packages("gratis")
```

### Development version

You can install the **development** version of `gratis` package from [GitHub
Repository](https://github.com/ykang/gratis) with:

``` r
devtools::install_github("ykang/gratis")
```

Usage
-----

### Load the package

``` r
require("gratis")
```

### Generate diverse time series

``` r
x <- generate_ts(n.ts = 2, freq = 12, nComp = 2, n = 120)
x$N1$pars
autoplot(x$N1$x)
```

### Generate mutiple seasonal time series

``` r
x <- generate_msts(seasonal.periods = c(7, 365), n = 800, nComp = 2)
autoplot(x)
```

### Generate time series with controllable features

``` r
x <- generate_ts_with_target(n = 1, ts.length = 60, freq = 1, seasonal = 0,
                             features = c('entropy', 'stl_features'),
                             selected.features = c('entropy', 'trend'),
                             target = c(0.6, 0.9))
autoplot(x)
```

### Web application

You could run the time series generation procedure in a web application
``` r
app_gratis()
```
Or visit our [online Shiny APP](https://ebsmonash.shinyapps.io/tsgeneration/)

See also
--------

- R package `tsfeatures` from [GitHub Repository](https://github.com/robjhyndman/tsfeatures).


References
----------

- Kang, Y., Hyndman, R., and Li, F. (2020). **GRATIS**: **G**ene**RA**ting **TI**me **S**eries with
diverse and controllable characteristics. [Statistical Analysis and Data Mining](https://doi.org/10.1002/sam.11461).


License
-------
This package is free and open source software, licensed under GPL-3.


Acknowledgements
----------------
Feng Li and Yanfei Kang are supported by the National Natural Science Foundation of China
(No. 11501587 and No. 11701022 respectively). Rob J Hyndman is supported by the Australian
Centre of Excellence in Mathematical and Statistical Frontiers.
