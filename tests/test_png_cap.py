"""P3: PNG output (the plot kit) is capped at 150 dpi and 3000 px per side."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from boldgenotyper import utils


def test_png_dpi_caps_resolution_and_size():
    fig = plt.figure(figsize=(6, 4))
    assert utils.png_dpi(fig, 300) == 150
    wide = plt.figure(figsize=(40, 5))
    assert utils.png_dpi(wide, 300) * (40 + 0.2) <= 3000
    plt.close("all")
