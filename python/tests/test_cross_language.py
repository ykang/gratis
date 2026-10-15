"""Compare fixed models against the repository's R implementation."""

import io
import itertools
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

import gratis
from gratis.target import _params_to_mar


ETS_CASES = [
    dict(error=error, seasonal=seasonal, trend=trend, sigma=0.05)
    for error, seasonal, trend in itertools.product(
        ["A", "M"], ["N", "A", "M"], ["N", "A", "Ad"]
    )
] + [dict(error="M", seasonal="N", trend="N", sigma=2.0)]

ARIMA_CASES = [
    dict(frequency=1, p=1, d=0, q=1, P=0, D=0, Q=0, phi=[0.99], theta=[0.3],
         Phi=[], Theta=[], constant=0, sigma=1e-8),
    dict(frequency=4, p=0, d=0, q=0, P=2, D=0, Q=0, phi=[], theta=[],
         Phi=[0.5, -0.5], Theta=[], constant=2, sigma=0.2),
    dict(frequency=4, p=0, d=0, q=0, P=0, D=0, Q=2, phi=[], theta=[],
         Phi=[], Theta=[0.5, -0.5], constant=-1, sigma=2),
    dict(frequency=1, p=6, d=0, q=0, P=0, D=0, Q=0,
         phi=[5.4, -12.15, 14.58, -9.8415, 3.54294, -.531441], theta=[],
         Phi=[], Theta=[], constant=0, sigma=1),
    dict(frequency=4, p=0, d=0, q=0, P=6, D=0, Q=0, phi=[], theta=[],
         Phi=[5.4, -12.15, 14.58, -9.8415, 3.54294, -.531441], Theta=[],
         constant=0, sigma=1),
] + [
    dict(frequency=4, p=1, d=d, q=1, P=2, D=D, Q=2, phi=[0.3], theta=[0.2],
         Phi=[0.5, -0.5], Theta=[0.2, -0.2], constant=2 if d + D <= 1 else 0,
         sigma=0.5)
    for d, D in [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
]

SCALE_CASES = [[1e155, 2e155, 3e155], [1000001, 1000002, 1000003],
               [-1e308, 0, 1e308], [1e-300, 2e-300, 3e-300], [5, 5, 5]]
WEIGHT_CASES = [[1], [.9, .05, .05], [.05, .9, .05], [0, 0, 1], [1, 0, 0], [.2, 0, .8]]


def _arima_innovations(case):
    generator = np.random.default_rng(7)
    return generator.normal(0, case["sigma"], 17), generator.normal(0, case["sigma"], 12)


def _r_vector(values):
    return "c(" + ",".join(format(float(x), ".17g") for x in values) + ")"


def _ets_parameters(case):
    return dict(
        frequency=4, error=case["error"], trend="N" if case["trend"] == "N" else "A",
        seasonal=case["seasonal"], damped=case["trend"] == "Ad",
        alpha=0.2, beta=0.05, gamma=0.1, phi=0.9, level=10, slope=0.2,
        season=[-1, 0.5, 1, -0.5] if case["seasonal"] == "A" else [0.7, 1.1, 1.3, 0.9],
        sigma=case["sigma"],
    )


@pytest.fixture(scope="module")
def r_results():
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript is needed for cross-language tests")
    script = [
        "if (!requireNamespace('forecast',quietly=TRUE) || "
        "!requireNamespace('polynom',quietly=TRUE)) quit(status=42)",
        "suppressPackageStartupMessages(library(forecast))",
        "library(polynom)",
        "source('R/arima_model.R'); source('R/simulate_arima.R'); "
        "source('R/ets_model.R'); source('R/mar_model.R'); source('R/simulate.R'); "
        "source('R/fitness_ts.R'); source('R/simulate_target.R'); source('R/rmixmorm.R')",
        "emit <- function(name, x) cat(name, ':', paste(sprintf('%.17g', x), "
        "collapse=','), '\\n', sep='')",
        "emit('pi_ma', pi_coefficients(ar=.3,ma=.5))",
        "emit('pi_seasonal', pi_coefficients(ar=.3,ma=.2,sar=.1,sma=.4,d=1,m=4))",
        "m <- mar_model(phi=matrix(.3),Phi=matrix(c(.2,.3)),d=0,D=0,constants=0,"
        "sigmas=1,weights=1,seasonal_periods=365); emit('mar_annual',m$ar)",
    ]
    for i, case in enumerate(ETS_CASES):
        parameters = _ets_parameters(case)
        innovations = np.random.default_rng(1).normal(0, case["sigma"], 12)
        script.append(
            f"m <- ets_model(frequency=4,error='{case['error']}',"
            f"trend='{parameters['trend']}',seasonal='{case['seasonal']}',"
            f"damped={'TRUE' if parameters['damped'] else 'FALSE'},alpha=.2,beta=.05,"
            f"gamma=.1,phi=.9,level=10,slope=.2,season={_r_vector(parameters['season'])},"
            f"sigma={case['sigma']}); emit('ets_{i}',stats::simulate(m,nsim=12,"
            f"innov={_r_vector(innovations)}))"
        )
    for i, (d, D) in enumerate([(0, 0), (1, 0), (0, 1)]):
        script.append(
            f"m <- arima_model(frequency=4,p=1,d={d},q=0,P=1,D={D},Q=0,"
            "phi=.5,Phi=.2,constant=2,sigma=1); "
            f"b <- arima_model(frequency=4,p=1,d={d},q=0,P=1,D={D},Q=0,"
            "phi=.5,Phi=.2,constant=0,sigma=1); "
            "set.seed(1); x <- stats::simulate(m,nsim=12,innov=rep(0,12)); "
            "set.seed(1); y <- stats::simulate(b,nsim=12,innov=rep(0,12)); "
            f"emit('arima_{i}',x-y)"
        )
    for i, case in enumerate(ARIMA_CASES):
        arguments = ",".join(f"{key}={_r_vector(value) if isinstance(value, list) else value}"
                             for key, value in case.items())
        startup, innovations = _arima_innovations(case)
        script.append(
            f"m <- arima_model({arguments}); "
            f"emit('arima_path_{i}',stats::simulate(m,innov={_r_vector(innovations)},"
            f"start.innov={_r_vector(startup)})); "
            f"emit('arima_zero_{i}',stats::simulate(m,innov={_r_vector(innovations)},n.start=0)); "
            f"emit('arima_burnin_{i}',arima_burnin(m))"
        )
    for i, values in enumerate(SCALE_CASES):
        script.append(f"emit('scaled_{i}',scalets({_r_vector(values)}))")
    for kind in ("arima", "ets"):
        for i, sigma in enumerate([1e155, 1e-200]):
            script.append(
                f"emit('sigma_{kind}_{i}',inherits(try({kind}_model(sigma={sigma}),"
                "silent=TRUE),'try-error'))"
            )
    script.extend([
        "original_mar_model <- mar_model",
        "mar_model <- function(...) { m <- original_mar_model(...); "
        "captured_weights <<- m$weights; m }",
    ])
    for i, weights in enumerate(WEIGHT_CASES):
        k = len(weights)
        parameters = np.r_[np.zeros(6 * k), np.ones(k), np.cumsum(weights)[:-1]]
        script.append(
            f"invisible(fitness_mar({_r_vector(parameters)},length=12,seasonal_periods=1,"
            f"ncomponents={k},feature_function=mean,target=0)); emit('target_weights_{i}',captured_weights)"
        )
    script.append("mar_model <- original_mar_model")
    script.extend([
        "emit('conditional_mixture',dmixnorm_ts(c(.2,-.3,.5,1,-2),"
        "list(c(.1,.2),c(-.1,.4,.1)),list(c(.5,1,1.5,2,2.5),rep(2,5)),c(2,3),log=TRUE))",
        "emit('conditional_inactive',dmixnorm_ts(c(1e308,1e308,0),"
        "list(c(0,0,0),c(0,2,-2)),list(rep(1,3),rep(1,3)),c(1,0),log=TRUE))",
        "emit('mixture_tail',dmixnorm(c(0,40,-3),c(0,1),c(1,2),c(2,3),log=TRUE))",
        "emit('mixture_zero_weight',dmixnorm(c(0,40),c(0,1),c(1,2),c(0,1),log=TRUE))",
        "emit('mixture_multivariate',dmixnorm(rbind(c(0,0),c(1,2)),matrix(c(.2,-.1),2,1),"
        "array(c(2,.3,.3,1),c(2,2,1)),2,log=TRUE))",
    ])
    root = Path(__file__).resolve().parents[2]
    result = subprocess.run(["Rscript", "--vanilla", "-"], input="\n".join(script),
                            cwd=root, capture_output=True, text=True, timeout=120)
    if result.returncode == 42:
        pytest.skip("R packages forecast and polynom are needed for cross-language tests")
    assert result.returncode == 0, result.stderr
    return {name: np.loadtxt(io.StringIO(values), delimiter=",", ndmin=1)
            for name, values in (line.split(":", 1) for line in result.stdout.splitlines())}


def test_ar_expansions_and_annual_mar_match_r(r_results):
    np.testing.assert_allclose(gratis.pi_coefficients(ar=0.3, ma=0.5), r_results["pi_ma"])
    np.testing.assert_allclose(
        gratis.pi_coefficients(ar=0.3, ma=0.2, sar=0.1, sma=0.4, d=1, m=4),
        r_results["pi_seasonal"],
    )
    model = gratis.mar_model(phi=[0.3], Phi=[0.2, 0.3], d=0, D=0, constants=[0],
                             sigmas=[1], weights=[1], seasonal_periods=365)
    np.testing.assert_allclose(model.ar[:, 0], r_results["mar_annual"])


@pytest.mark.parametrize("index,case", list(enumerate(ETS_CASES)))
def test_all_ets_components_match_r_under_identical_innovations(r_results, index, case):
    model = gratis.ets_model(**_ets_parameters(case))
    np.testing.assert_allclose(model.simulate(12, rng=1), r_results[f"ets_{index}"],
                               rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("backend", ["statsmodels", "numpy"])
@pytest.mark.parametrize("index,orders", list(enumerate([(0, 0), (1, 0), (0, 1)])))
def test_arima_mean_and_drift_effects_match_r(r_results, backend, index, orders):
    d, D = orders
    parameters = dict(frequency=4, p=1, d=d, q=0, P=1, D=D, Q=0,
                      phi=[0.5], Phi=[0.2], sigma=1, backend=backend)
    model = gratis.arima_model(**parameters, constant=2)
    baseline = gratis.arima_model(**parameters, constant=0)
    difference = model.simulate(12, rng=1) - baseline.simulate(12, rng=1)
    np.testing.assert_allclose(difference, r_results[f"arima_{index}"], atol=1e-12)


@pytest.mark.parametrize("backend", ["statsmodels", "numpy"])
@pytest.mark.parametrize("index,case", list(enumerate(ARIMA_CASES)))
def test_arima_complete_paths_and_burnin_match_r(r_results, backend, index, case):
    model = gratis.arima_model(**case, backend=backend)
    startup, innovations = _arima_innovations(case)
    np.testing.assert_allclose(model.simulate(innov=innovations, start_innov=startup),
                               r_results[f"arima_path_{index}"], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(model.simulate(innov=innovations, n_start=0),
                               r_results[f"arima_zero_{index}"], rtol=1e-12, atol=1e-12)
    assert model._default_burnin() == r_results[f"arima_burnin_{index}"][0]


@pytest.mark.parametrize("kind", ["arima", "ets"])
@pytest.mark.parametrize("index,sigma", list(enumerate([1e155, 1e-200])))
def test_noise_variance_validation_matches_r(r_results, kind, index, sigma):
    assert r_results[f"sigma_{kind}_{index}"][0] == 1
    with pytest.raises(ValueError, match="normal variance"):
        getattr(gratis, f"{kind}_model")(sigma=sigma)


@pytest.mark.parametrize("index,values", list(enumerate(SCALE_CASES)))
def test_scaling_matches_r_for_extreme_and_offset_inputs(r_results, index, values):
    np.testing.assert_allclose(gratis.scalets(values), r_results[f"scaled_{index}"], atol=1e-12)


@pytest.mark.parametrize("index,weights", list(enumerate(WEIGHT_CASES)))
def test_target_weight_parameterization_matches_r(r_results, index, weights):
    k = len(weights)
    parameters = np.r_[np.zeros(6 * k), np.ones(k), np.cumsum(weights)[:-1]]
    model = _params_to_mar(parameters, np.array([1]), k)
    np.testing.assert_allclose(model.weights, r_results[f"target_weights_{index}"], atol=1e-15)


def test_mixture_densities_match_r_after_normalization_and_in_tails(r_results):
    np.testing.assert_allclose(gratis.dmixnorm([0, 40, -3], [0, 1], [1, 2], [2, 3], log=True),
                               r_results["mixture_tail"], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(gratis.dmixnorm([0, 40], [0, 1], [1, 2], [0, 1], log=True),
                               r_results["mixture_zero_weight"], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(
        gratis.dmixnorm([[0, 0], [1, 2]], np.array([[.2], [-.1]]),
                        np.array([[2, .3], [.3, 1]])[:, :, None], [2], log=True),
        r_results["mixture_multivariate"], rtol=1e-12, atol=1e-12)


def test_vectorized_conditional_densities_match_r_including_inactive_components(r_results):
    np.testing.assert_allclose(
        gratis.dmixnorm_ts([.2, -.3, .5, 1, -2], [[.1, .2], [-.1, .4, .1]],
                           [[.5, 1, 1.5, 2, 2.5], [2]*5], [2, 3], log=True),
        r_results["conditional_mixture"], rtol=1e-12, atol=1e-12)
    with np.errstate(all="raise"):
        actual = gratis.dmixnorm_ts([1e308, 1e308, 0], [[0, 0, 0], [0, 2, -2]],
                                    [[1]*3, [1]*3], [1, 0], log=True)
    np.testing.assert_allclose(actual, r_results["conditional_inactive"], atol=1e-12)
