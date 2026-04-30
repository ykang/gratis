import numpy as np
import pytest

import gratis


class FakeAxes:
    def __init__(self):
        self.lines = []
        self.title = None
        self.xlabel = None
        self.ylabel = None
        self.legend_called = False

    def plot(self, x, y, **kwargs):
        self.lines.append((np.asarray(x), np.asarray(y), kwargs))

    def set_title(self, title):
        self.title = title

    def set_xlabel(self, label):
        self.xlabel = label

    def set_ylabel(self, label):
        self.ylabel = label

    def legend(self):
        self.legend_called = True


def test_plot_series_plots_single_series_on_supplied_axes():
    ax = FakeAxes()

    returned = gratis.plot_series(np.array([1.0, 2.0, 3.0]), ax=ax, title="Generated")

    assert returned is ax
    assert len(ax.lines) == 1
    np.testing.assert_array_equal(ax.lines[0][0], np.array([0, 1, 2]))
    np.testing.assert_array_equal(ax.lines[0][1], np.array([1.0, 2.0, 3.0]))
    assert ax.lines[0][2]["label"] == "Series 1"
    assert ax.title == "Generated"
    assert ax.xlabel == "Index"
    assert ax.ylabel == "Value"
    assert not ax.legend_called


def test_plot_series_plots_multiple_series_with_labels():
    ax = FakeAxes()
    data = np.column_stack([np.arange(4), np.arange(4) + 10])

    gratis.plot_series(data, ax=ax, index=[10, 11, 12, 13], labels=["a", "b"])

    assert len(ax.lines) == 2
    np.testing.assert_array_equal(ax.lines[0][0], np.array([10, 11, 12, 13]))
    assert ax.lines[0][2]["label"] == "a"
    assert ax.lines[1][2]["label"] == "b"
    assert ax.legend_called


def test_plot_series_respects_max_series():
    ax = FakeAxes()
    data = np.ones((5, 3))

    gratis.plot_series(data, ax=ax, max_series=2)

    assert len(ax.lines) == 2


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"data": np.ones((2, 2, 2))}, "one-dimensional series"),
        ({"data": np.ones((2, 2)), "labels": ["one"]}, "labels"),
        ({"data": np.ones(2), "index": [1]}, "index length"),
        ({"data": np.ones((2, 2)), "max_series": 0}, "max_series"),
    ],
)
def test_plot_series_validates_inputs(kwargs, message):
    with pytest.raises(ValueError, match=message):
        gratis.plot_series(ax=FakeAxes(), **kwargs)


def test_plot_generated_generates_and_plots_model():
    ax = FakeAxes()
    model = gratis.arima_model(p=1, d=0, q=0, phi=[0.2], constant=0, sigma=1)

    returned = gratis.plot_generated(model, length=8, nseries=2, ax=ax, rng=1)

    assert returned is ax
    assert len(ax.lines) == 2
    assert ax.lines[0][1].shape == (8,)
