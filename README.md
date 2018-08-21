`tsgeneration`
==============

The R package `tsgeneration` provides efficient algorithms for generating time series with
diverse and controllable characteristics.

Installation
------------

The package requires the [`tsfeatures`](https://github.com/robjhyndman/tsfeatures) package
to be installed.

``` r
# install.packages("devtools")
devtools::install_github("robjhyndman/tsfeatures")
```

You can install the **development** version of `tsgeneration` package from [GitHub
Repository](https://github.com/ykang/tsgeneration) with:

``` r
devtools::install_github("ykang/tsgeneration")
```

Usage
-----

### Load the package

``` r
require("tsgeneration")
```

### Generate diverse time series

``` r
x <- generate_ts(n.ts = 2, freq = 12, nComp = 2, n = 120)
x$N1$pars
plot(x$N1$x)
```

### Simulate mutiple seasonal time series

``` r
x <- simulate_msts(seasonal.periods = c(7, 365), n = 800, nComp = 2)
plot(x)
```

### Web application

You could run the time series generation procedure in a web application
``` r
run_tsgeneration_app()
```

See also
--------

- R package ['tsfeatures'](https://github.com/robjhyndman/tsfeatures)


References
----------

- Kang, Y., Hyndman, R., and Li, F. (2018) _Efficient generation of time series with
diverse and controllable characteristics_ Working paper.


License
-------
This package is free and open source software, licensed under GPL-3.
