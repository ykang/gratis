`tsgeneration`
==============

The R package `tsgeneration` provides efficient algorithms for generating time series with
diverse and controllable characteristics.

Installation
------------

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
run_tsgeneration_app()
```

See also
--------

- R package `tsfeatures` from [GitHub Repository](https://github.com/robjhyndman/tsfeatures).


References
----------

- Kang, Y., Hyndman, R., and Li, F. (2018). Efficient generation of time series with
diverse and controllable characteristics. [Working paper](https://robjhyndman.com/publications/tsgeneration/).


License
-------
This package is free and open source software, licensed under GPL-3.


Acknowledgements
----------------
Feng Li and Yanfei Kang's research were supported by the
National Natural Science Foundation of China
(No. 11501587 and No. 11701022 respectively).
