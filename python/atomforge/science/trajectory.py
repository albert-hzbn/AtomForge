"""Multi-frame import/export and GIF rendering through ASE and Matplotlib."""

from pathlib import Path
from .simulation import to_ase, from_ase


def read_trajectory(path, *, format=None, index=":"):
    """Read XYZ/extXYZ, XDATCAR, LAMMPS dumps and other ASE-supported trajectories."""
    from ase.io import read
    atoms = read(str(path), index=index, format=format)
    return [from_ase(frame) for frame in (atoms if isinstance(atoms,list) else [atoms])]


def write_trajectory(frames, path, *, format=None):
    """Write a nonempty sequence of structures using an ASE trajectory format."""
    from ase.io import write
    frames = list(frames)
    if not frames:
        raise ValueError("Trajectory is empty")
    write(str(path), [to_ase(frame) for frame in frames], format=format)


def export_animation(frames, path, *, fps=12, rotation="20x,30y", radii=0.5):
    """Export a GIF with fixed spatial bounds and explicitly chosen atom radii."""
    import numpy as np
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.animation import FuncAnimation, PillowWriter
    from ase.visualize.plot import plot_atoms
    frames = [to_ase(frame) for frame in frames]
    if not frames or fps <= 0 or not np.isfinite(fps) or radii <= 0:
        raise ValueError("Provide frames and positive fps/radii")
    if Path(path).suffix.lower() != ".gif":
        raise ValueError("Animation output must use .gif")
    figure = Figure(figsize=(7,7)); FigureCanvasAgg(figure)
    axes = figure.add_subplot(111)
    # Determine bounds in the requested projection across every frame.
    bounds = []
    for frame in frames:
        axes.clear(); plot_atoms(frame, axes, rotation=rotation, radii=radii)
        bounds.append((*axes.get_xlim(), *axes.get_ylim()))
    bounds = np.asarray(bounds)
    xlim = (bounds[:,0].min(), bounds[:,1].max())
    ylim = (bounds[:,2].min(), bounds[:,3].max())
    def draw(index):
        axes.clear(); plot_atoms(frames[index], axes, rotation=rotation, radii=radii)
        axes.set_xlim(xlim); axes.set_ylim(ylim); axes.set_axis_off()
        axes.set_title("Frame {} / {}".format(index+1,len(frames)))
    animation = FuncAnimation(figure, draw, frames=len(frames), interval=1000/fps)
    animation.save(str(path), writer=PillowWriter(fps=fps))
