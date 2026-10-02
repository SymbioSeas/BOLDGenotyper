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


def test_point_sizes_when_every_location_has_the_same_count():
    # min-max scaling divided by zero here: NaN sizes broke the legend at save
    import pandas as pd
    sizes = visualization._point_sizes(pd.Series([2, 2, 2]), 50, 500)
    assert sizes.notna().all() and (sizes == 50).all()
    sizes = visualization._point_sizes(pd.Series([1, 3, 5]), 50, 500)
    assert sizes.tolist() == [50.0, 275.0, 500.0]


def test_map_saves_when_every_location_has_the_same_count(tmp_path, caplog):
    import logging
    import pandas as pd
    # three locations, two records each (counts all 2)
    df = pd.DataFrame({
        "haplotype": ["H1", "H1", "H2", "H2", "H1", "H1"],
        "lat": [34.36, 34.36, 38.98, 38.98, 38.98, 38.98],
        "lon": [-120.0, -120.0, -123.95, -123.95, -123.95, -123.95],
    })
    out = tmp_path / "map.pdf"
    with caplog.at_level(logging.WARNING, logger="boldgenotyper"):
        visualization.plot_distribution_map(df=df, output_path=str(out), genotype_column="haplotype",
                                            latitude_col="lat", longitude_col="lon")
    assert out.exists() and out.read_bytes().startswith(b"%PDF")
    assert "Failed to save" not in caplog.text
