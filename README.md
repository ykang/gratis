`tsgeneration`
==============

The R package `tsgeneration` provides efficient algorithms for generating time series with
diverse and controllable characteristics.

Installation
------------

You can install the **development** version from [GitHub
Repository](https://github.com/ykang/tsgeneration) with:

``` r
# install.packages("devtools")
devtools::install_github("ykang/tsgeneration")

```

Usage
-----

### Load the package

``` r
require("tsgeneration")
```


### Web Application

You could run the time series generation procedure in a web application
``` r
tsgeneration:::run.tsgeneration.app()
```


References
----------

- Kang, Y., Hyndman, R., and Li, F. (2018) _Efficient generation of time series with
diverse and controllable characteristics_ Working paper.


License
-------
This package is free and open source software, licensed under GPL-3.
