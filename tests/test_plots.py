"""Test integrated plots."""
from pathlib import Path

import numpy as np

from taurenmd.plots import labeldots, plotparams


def test_plotparams_0():
    """Test plotparams defaults."""
    plotparams.plot(
        list(range(100)),
        list(range(100)),
        )

    Path('plot_param.pdf').unlink()
    return


def test_plotparams_1():
    """Test plotparams swapped arguments."""
    plotparams.plot(
        list(range(100)),
        list(range(100)),
        labels='This labels',
        title='Title',
        xlabel='data',
        ylabel='observable',
        colors=['green'],
        xmax=50,
        xmin=20,
        ymax=50,
        ymin=20,
        grid=False,
        legend=False,
        filename='this_plot.png',
        )

    Path('this_plot.png').unlink(missing_ok=False)
    return


def test_plotparams_2():
    """Test plotparams swapped arguments with multidimensional y-data."""
    plotparams.plot(
        list(range(50)),
        np.arange(100).reshape(2, 50),
        labels=['label 1', 'label 2'],
        title='Title',
        xlabel='data',
        ylabel='observable',
        colors=['green', 'blue'],
        xmax=50,
        xmin=20,
        ymax=50,
        ymin=20,
        grid=False,
        legend=False,
        vert_lines=(1, 10),
        filename='this_plot.png',
        )

    Path('this_plot.png').unlink(missing_ok=False)
    return


def test_labeldots_0():
    """Test label dots defaults."""
    labeldots.plot(
        list(range(50)),
        )

    Path('plot_param.pdf').unlink(missing_ok=False)


def test_labeldots_1():
    """Test label dots defaults."""
    labeldots.plot(
        list(range(50)),
        x_labels=[f'{i}_' for i in range(50)],
        xlabel='x label',
        x_label_rot=45,
        labels=['series A'],
        title='title',
        legend=False,
        grid=False,
        filename='plot_test.png',
        )

    Path('plot_test.png').unlink(missing_ok=False)


def test_labeldots_2():
    """Test label dots defaults."""
    labeldots.plot(
        np.arange(100).reshape(2, 50),
        x_labels=[f'{i}_' for i in range(50)],
        x_label_rot=45,
        colors=['k', 'w'],
        labels=['series A', 'series B'],
        title='title',
        legend=False,
        grid=False,
        filename='plot_test.png',
        )

    Path('plot_test.png').unlink(missing_ok=False)
