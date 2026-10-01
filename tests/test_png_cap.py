"""P3: PNG output is capped at 150 dpi and 3000 px per side; metadata figures are PDF + PNG."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

from boldgenotyper import metadata_summary as ms
from boldgenotyper import utils


def test_png_dpi_caps_resolution_and_size():
    fig = plt.figure(figsize=(6, 4))
    assert utils.png_dpi(fig, 300) == 150
    wide = plt.figure(figsize=(40, 5))
    assert utils.png_dpi(wide, 300) * (40 + 0.2) <= 3000
    plt.close("all")


def test_metadata_figures_are_pdf_and_capped_png(tmp_path):
    out = tmp_path / "metadata" / "pdf" / "wide"
    out.parent.mkdir(parents=True)
    fig = plt.figure(figsize=(40, 5))
    plt.plot([0, 1], [0, 1])
    ms._save_multi_format(fig, out, dpi=300)
    plt.close("all")
    assert (tmp_path / "metadata" / "pdf" / "wide.pdf").exists()
    assert not (tmp_path / "metadata" / "svg").exists()
    png = Image.open(tmp_path / "metadata" / "png" / "wide.png")
    assert max(png.size) <= 3000
    assert png.info.get("dpi", (150,))[0] <= 150.5
