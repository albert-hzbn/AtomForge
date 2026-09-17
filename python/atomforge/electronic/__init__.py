"""Electronic post-processing with AtomForge's native calculation library.

The library is loaded lazily; importing this module needs no native binaries.
"""

from ._grid import Grid, Surface, Volume, BaderPartition, load_volume
from ._models import miller_section, model_density, reciprocal_vectors, scattering_factor
from .pipeline import apply_operation, run_pipeline, batch_process
from ._wannier import WannierHamiltonian, read_hr
from ._lobster import read_cohpcar, read_icohplist

__all__ = ["Grid", "Surface", "Volume", "BaderPartition", "load_volume", "miller_section",
           "model_density", "reciprocal_vectors", "scattering_factor",
           "apply_operation", "run_pipeline", "batch_process",
           "WannierHamiltonian", "read_hr", "read_cohpcar", "read_icohplist"]
