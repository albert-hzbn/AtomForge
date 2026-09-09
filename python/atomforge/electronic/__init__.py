"""Electronic post-processing with AtomForge's native calculation library.

The library is loaded lazily; importing this module needs no native binaries.
See ELECTRONIC_POSTPROCESSING.md for build instructions, units and methods.
"""

from ._grid import Grid, Surface, Volume, load_volume
from ._models import miller_section, model_density, reciprocal_vectors, scattering_factor

__all__ = ["Grid", "Surface", "Volume", "load_volume", "miller_section",
           "model_density", "reciprocal_vectors", "scattering_factor"]
