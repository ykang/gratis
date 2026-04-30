import numpy as np
import pytest

import gratis


def _require_statsmodels():
    try:
        return pytest.importorskip("statsmodels")
    except Exception as exc:
        pytest.skip(f"statsmodels is not importable: {exc}")


def test_public_imports():
    assert gratis.__version__ == "0.0.0"
    assert callable(gratis.rmixnorm)
    assert callable(gratis.dmixnorm)
    assert callable(gratis.pi_coefficients)
    assert callable(gratis.arima_model)
    assert callable(gratis.mar_model)
    assert callable(gratis.ets_model)


def test_rmixnorm_supports_more_than_five_components():
    draws = gratis.rmixnorm(
        100,
        means=np.arange(6, dtype=float),
        sigmas=np.ones(6),
        weights=np.full(6, 1 / 6),
        rng=1,
    )

    assert draws.shape == (100,)
    assert np.isfinite(draws).all()


def test_pi_coefficients_basic_cases():
    np.testing.assert_allclose(gratis.pi_coefficients(ar=0.8), np.array([0.8]))
    np.testing.assert_allclose(
        gratis.pi_coefficients(ar=[0.4, -0.2], d=1),
        np.array([1.4, -0.6, 0.2]),
    )


def test_dmixnorm_returns_density_and_log_density():
    density = gratis.dmixnorm(
        np.array([-1.0, 0.0, 1.0]),
        means=np.array([0.0]),
        sigmas=np.array([1.0]),
        weights=np.array([1.0]),
    )
    log_density = gratis.dmixnorm(
        np.array([-1.0, 0.0, 1.0]),
        means=np.array([0.0]),
        sigmas=np.array([1.0]),
        weights=np.array([1.0]),
        log=True,
    )

    np.testing.assert_allclose(np.log(density), log_density)


def test_mar_model_constructs_and_simulates_explicit_model():
    phi = np.array([[0.8, 0.6], [0.0, 0.3]])
    model = gratis.mar_model(
        phi=phi,
        d=0,
        sigmas=np.array([1.0, 2.0]),
        weights=np.array([0.8, 0.2]),
    )

    assert model.k == 2
    assert model.ar.shape == (3, 2)
    draws = model.simulate(n=50, rng=1)
    assert draws.shape == (50,)
    assert np.isfinite(draws).all()


def test_mar_model_handles_zero_order_zero_constant_model():
    model = gratis.mar_model(
        k=1,
        p=0,
        d=0,
        constants=0,
        sigmas=1,
        weights=[1],
    )

    assert model.ar.shape == (1, 1)
    assert model.ar[0, 0] == 0


def test_arima_model_constructs_and_simulates():
    model = gratis.arima_model(
        p=1,
        d=0,
        q=1,
        phi=[0.4],
        theta=[0.2],
        constant=0,
        sigma=1,
    )

    assert model.order == (1, 0, 1)
    draws = gratis.simulate(model, n=40, rng=1)
    assert draws.shape == (40,)
    assert np.isfinite(draws).all()


def test_arima_model_can_require_statsmodels_backend():
    _require_statsmodels()
    model = gratis.arima_model(
        p=1,
        d=0,
        q=1,
        phi=[0.4],
        theta=[0.2],
        constant=0,
        sigma=1,
        backend="statsmodels",
    )

    draws = model.simulate(n=20, rng=1)
    assert draws.shape == (20,)
    assert np.isfinite(draws).all()


def test_ets_model_constructs_and_simulates():
    model = gratis.ets_model(
        error="A",
        trend="N",
        seasonal="N",
        alpha=0.2,
        level=1.0,
        sigma=1.0,
    )

    assert model.method == "ETS(A,N,N)"
    draws = model.simulate(n=40, rng=1)
    assert draws.shape == (40,)
    assert np.isfinite(draws).all()


def test_ets_model_can_require_statsmodels_backend():
    _require_statsmodels()
    model = gratis.ets_model(
        error="A",
        trend="N",
        seasonal="N",
        alpha=0.2,
        level=1.0,
        sigma=1.0,
        backend="statsmodels",
    )

    draws = model.simulate(n=20, rng=1)
    assert draws.shape == (20,)
    assert np.isfinite(draws).all()


def test_generate_returns_distinct_series_when_seeded_once():
    model = gratis.arima_model(p=1, d=0, q=0, phi=[0.2], constant=0, sigma=1)
    generated = gratis.generate(model, length=20, nseries=2, rng=1)

    assert generated.shape == (20, 2)
    assert not np.array_equal(generated[:, 0], generated[:, 1])


def test_simulate_target_runs_lightweight_search():
    series = gratis.simulate_target(
        length=12,
        feature_function=lambda x: np.array([np.mean(x)]),
        target=np.array([0.0]),
        tolerance=10.0,
        population_size=4,
        max_iter=1,
        rng=1,
    )

    assert series.shape == (12,)
    assert np.isfinite(series).all()
