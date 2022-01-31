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

    assert Path('plot_param.pdf').exists()
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

    # py37 compatibility
    assert Path('this_plot.png').exists()
    Path('this_plot.png').unlink()
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

    # py37 compatibility
    assert Path('this_plot.png').exists()
    Path('this_plot.png').unlink()
    return


def test_plotparams_3():
    """Test plotparams swapped arguments with multidimensional y-data."""
    plotparams.plot(
        list(range(50)),
        np.arange(100).reshape(2, 50),
        labels=['label 1', 'label 2'],
        title='Title',
        xlabel='data',
        xticks=list(range(50)),
        xticks_labels=list(map(str, range(50))),
        yticks=list(range(10)),
        yticks_labels=list(map(str, range(10))),
        hline=0.5,
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

    # py37 compatibility
    assert Path('this_plot.png').exists()
    Path('this_plot.png').unlink()
    return


def test_labeldots_0():
    """Test label dots defaults."""
    labeldots.plot(
        list(range(50)),
        )

    # py37 compatibility
    assert Path('plot_param.pdf').exists()
    Path('plot_param.pdf').unlink()


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

    # py37 compatibility
    assert Path('plot_test.png').exists()
    Path('plot_test.png').unlink()


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

    # py37 compatibility
    assert Path('plot_test.png').exists()
    Path('plot_test.png').unlink()
