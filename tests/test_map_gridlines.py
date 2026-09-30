"""P0-6: the main distribution map failed with 'dash list must be positive'."""

import io

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

ccrs = pytest.importorskip("cartopy.crs")

from boldgenotyper import visualization


def test_main_map_gridlines_render_to_pdf():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_global()
    visualization._add_map_gridlines(ax)
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf")
    plt.close(fig)
    assert buf.getvalue().startswith(b"%PDF")
