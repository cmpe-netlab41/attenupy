"""Geometry builders: produce the PENGEOM ``.geo`` file for a simulation.

A geometry builder is any object with::

    filename(spec) -> str          # e.g. "oneL015.geo" (<= 80 chars)
    write(spec, path) -> None      # write the file for this spec
    validate(spec) -> None         # optional: raise ValueError if unsupported

Material numbering is a contract shared with the penEasy deck
(:mod:`attenupy.template`): 1 = shielding slab, 2 = air detector,
3 = air enclosure.

:class:`SlabGeometry` generates the default setup (a 9x9 cm slab at 100 cm
from a point source, a 4x4x1 cm air detector at 101 cm, all inside a 2 m air
sphere); every dimension is a parameter. :class:`FileGeometry` copies
user-provided files instead.
"""

import math
import os
import shutil
from dataclasses import dataclass
from decimal import Decimal

from .naming import mm_to_token

_SEP = "0" * 64


def _e(x):
    """PENGEOM value field, standard style: `` 1.000150000000000E+02``."""
    return f"{x: .15E}"


def _fortran_e(x):
    """PENGEOM value field, 0.xxx style: `` 0.450000000000000E+01``."""
    if x == 0:
        return " 0.000000000000000E+00"
    d = Decimal(repr(abs(float(x))))
    exp = d.adjusted() + 1
    return f"{'-' if x < 0 else ' '}{d.scaleb(-exp):.15f}E{exp:+03d}"


def _g(x):
    return f"{x:g}"


_SLAB_TEMPLATE = """\
{title_mm} mm shielding material
{width}x{width} cm2
{sep}
SURFACE (   1)   Z={z_front} cm
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=({z_front_e},   0)
{sep}
SURFACE (   2)   Z={z_back:.3f} cm ({thickness_mm:.2f} mm thick layer)
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=({z_back_e},   0)
{sep}
SURFACE (   3)   Y=-{hw} cm
INDICES=( 0, 0, 0, 0, 0)
     AY=( 1.000000000000000E+00,   0)
     A0=({hw_pos},   0)
{sep}
SURFACE (   4)   Y=+{hw} cm
INDICES=( 0, 0, 0, 0, 0)
     AY=( 1.000000000000000E+00,   0)
     A0=({hw_neg},   0)
{sep}
SURFACE (   5)   X=-{hw} cm
INDICES=( 0, 0, 0, 0, 0)
     AX=( 1.000000000000000E+00,   0)
     A0=({hw_pos},   0)
{sep}
SURFACE (   6)   X={hw} cm
INDICES=( 0, 0, 0, 0, 0)
     AX=( 1.000000000000000E+00,   0)
     A0=({hw_neg},   0)
{sep}
BODY    (   1)   Layer 1 (material)
MATERIAL(   1)
SURFACE (   1), SIDE POINTER=( 1)
SURFACE (   2), SIDE POINTER=(-1)
SURFACE (   3), SIDE POINTER=( 1)
SURFACE (   4), SIDE POINTER=(-1)
SURFACE (   5), SIDE POINTER=( 1)
SURFACE (   6), SIDE POINTER=(-1)
c
c
{sep}
SURFACE (  13)   Z = {det_front} cm
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=({det_front_e},   0)
{sep}
SURFACE (  14)   Z = {det_back} cm (detector {det_depth} cm thick)
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=({det_back_e},   0)
{sep}
SURFACE (  15)   X=-{dhw} cm
INDICES=( 0, 0, 0, 0, 0)
     AX=( 1.000000000000000E+00,   0)
     A0=({dhw_pos},   0)
{sep}
SURFACE (  16)   X={dhw} cm
INDICES=( 0, 0, 0, 0, 0)
     AX=( 1.000000000000000E+00,   0)
     A0=({dhw_neg},   0)
{sep}
SURFACE (  17)   Y=-{dhw} cm
INDICES=( 0, 0, 0, 0, 0)
     AY=( 1.000000000000000E+00,   0)
     A0=({dhw_pos},   0)
{sep}
SURFACE (  18)   Y={dhw} cm
INDICES=( 0, 0, 0, 0, 0)
     AY=( 1.000000000000000E+00,   0)
     A0=({dhw_neg},   0)
{sep}
BODY    (   2)  Air Detector, at {det_front} cm from the source, {gap} cm from the textile shield
MATERIAL(   2)
SURFACE (  13), SIDE POINTER=( 1)
SURFACE (  14), SIDE POINTER=(-1)
SURFACE (  15), SIDE POINTER=( 1)
SURFACE (  16), SIDE POINTER=(-1)
SURFACE (  17), SIDE POINTER=( 1)
SURFACE (  18), SIDE POINTER=(-1)
c
c
{sep}
SURFACE (  39)   sphere diameter={diameter}
INDICES=( 1, 1, 1, 0,-1)
X-SCALE=({radius_e},   0)
Y-SCALE=({radius_e},   0)
Z-SCALE=({radius_e},   0)
{sep}
BODY    (   3)  enclosure (air)
MATERIAL(   3)
BODY    (   1)
BODY    (   2)
SURFACE (  39), SIDE POINTER=(-1)
{sep}
END      0000000000000000000000000000000000000000000000000000000
"""


@dataclass(frozen=True)
class SlabGeometry:
    """A square slab facing a point source, followed by an air detector.

    All lengths in cm (except the slab thickness, which comes from the spec in
    mm). The source sits at the origin shooting along +z.
    """

    slab_z_cm: float = 100.0            # source -> slab front face
    slab_half_width_cm: float = 4.5     # slab is 9x9 cm
    detector_z_cm: float = 101.0        # source -> detector front face
    detector_depth_cm: float = 1.0
    detector_half_width_cm: float = 2.0  # detector is 4x4 cm
    world_radius_cm: float = 200.0      # air sphere enclosing everything

    @property
    def max_thickness_mm(self):
        return (self.detector_z_cm - self.slab_z_cm) * 10.0

    def validate(self, spec):
        if spec.thickness_mm >= self.max_thickness_mm:
            raise ValueError(
                f"Slab of {spec.thickness_mm} mm would overlap the detector; this "
                f"geometry allows < {self.max_thickness_mm:g} mm "
                "(increase detector_z_cm)"
            )

    def check_beam(self, half_angle_deg):
        """Raise if a source cone of this half-aperture is wider than the slab
        at its front face, i.e. part of the beam would bypass the shield."""
        radius = self.slab_z_cm * math.tan(math.radians(half_angle_deg))
        if radius > self.slab_half_width_cm + 1e-9:
            raise ValueError(
                f"The {half_angle_deg:g} deg source cone reaches {radius:.3f} cm from the axis "
                f"at the slab, but the slab half-width is {self.slab_half_width_cm:g} cm: part of the "
                "beam would reach the detector unshielded. Enlarge slab_half_width_cm or "
                "lower PeneasyDeck(beam_half_angle_deg=...)"
            )

    def filename(self, spec):
        try:
            return f"{mm_to_token(spec.thickness_mm)}.geo"
        except ValueError:  # not a multiple of 0.01 mm
            return f"slab_{spec.thickness_mm:g}mm.geo"

    def render(self, spec):
        self.validate(spec)
        mm = spec.thickness_mm
        z_back = round(self.slab_z_cm + mm / 10.0, 10)
        det_back = self.detector_z_cm + self.detector_depth_cm
        hw, dhw = self.slab_half_width_cm, self.detector_half_width_cm
        return _SLAB_TEMPLATE.format(
            sep=_SEP,
            title_mm=repr(float(mm)),
            width=_g(2 * hw),
            thickness_mm=mm,
            z_front=_g(self.slab_z_cm), z_front_e=_e(self.slab_z_cm),
            z_back=z_back, z_back_e=_e(z_back),
            hw=_g(hw), hw_pos=_fortran_e(hw), hw_neg=_fortran_e(-hw),
            det_front=_g(self.detector_z_cm), det_front_e=_e(self.detector_z_cm),
            det_back=_g(det_back), det_back_e=_e(det_back),
            det_depth=_g(self.detector_depth_cm),
            dhw=_g(dhw), dhw_pos=_fortran_e(dhw), dhw_neg=_fortran_e(-dhw),
            gap=_g(self.detector_z_cm - self.slab_z_cm),
            diameter=_g(2 * self.world_radius_cm),
            radius_e=_e(self.world_radius_cm),
        )

    def write(self, spec, path):
        with open(path, "w") as f:
            f.write(self.render(spec))


class FileGeometry:
    """Copy pre-made ``.geo`` files from ``directory``.

    ``pattern`` is formatted with ``token`` (``oneL015``) and ``thickness_mm``.
    """

    def __init__(self, directory, pattern="{token}.geo"):
        self.directory = os.path.abspath(directory)
        self.pattern = pattern

    def filename(self, spec):
        return self.pattern.format(token=mm_to_token(spec.thickness_mm),
                                   thickness_mm=spec.thickness_mm)

    def _source(self, spec):
        return os.path.join(self.directory, self.filename(spec))

    def validate(self, spec):
        if not os.path.isfile(self._source(spec)):
            raise FileNotFoundError(f"Geometry file not found: {self._source(spec)}")

    def write(self, spec, path):
        shutil.copyfile(self._source(spec), path)
