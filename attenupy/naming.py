"""Folder naming: how a :class:`SimulationSpec` becomes a directory name.

Since every simulation directory stores its exact parameters in ``spec.json``,
the folder name is only a human-readable, unique label. Any object with a
``name(spec) -> str`` method can be used as a namer (see :class:`Namer`).

The default :class:`CompactNamer` produces names like::

    RQR6_080_2_04_oneL080_BaSO4_Bi2O3_0.20_0.40
    |______| | |  |_____| |_________| |_______|
     beam    | |  thickness names       fractions (2 decimals)
            n  polymer*10 (oneL + mm*100)

and can parse them back, so a folder can be recognised even without its
``spec.json``.
"""

import re
from typing import Protocol

from .spec import Mixture, SimulationSpec

_FOLDER_RE = re.compile(r"^(.+?)_(\d+)_(\d{2})_oneL(\d+)_(.+)$")
_TOKEN_RE = re.compile(r"^oneL(\d+)$")


class Namer(Protocol):
    def name(self, spec):
        """Return the directory name (no path separators) for ``spec``."""


# ---------------------------------------------------------------------------
# Thickness tokens
# ---------------------------------------------------------------------------
def mm_to_token(thickness_mm):
    """0.15 -> ``"oneL015"``; 2.5 -> ``"oneL250"`` (hundredths of a mm)."""
    hundredths = round(thickness_mm * 100)
    if abs(hundredths - thickness_mm * 100) > 1e-6:
        raise ValueError(
            f"Thickness {thickness_mm} mm is not a multiple of 0.01 mm and has no "
            "'oneL' token"
        )
    return f"oneL{hundredths:03d}"


def token_to_mm(token):
    """``"oneL015"`` -> 0.15."""
    m = _TOKEN_RE.match(token)
    if not m:
        raise ValueError(f"{token!r} is not a thickness token like 'oneL015'")
    return int(m.group(1)) / 100.0


def parse_folder_name(name):
    """Decode a :class:`CompactNamer` folder name without building a spec.

    Returns a dict with keys: ``energy``, ``n_materials``, ``poly_ratio``,
    ``thickness`` (the ``oneL...`` token), ``thickness_mm`` (int, in
    hundredths of a mm -- the ``geo_number_mm`` column), ``materials``.

    Raises :class:`ValueError` if the name does not match the convention.
    """
    m = _FOLDER_RE.match(name)
    if not m:
        raise ValueError(f"{name!r} does not match the folder convention")
    energy, n_str, poly_str, geo_str, rest = m.groups()
    n_mat = int(n_str)
    parts = rest.split("_")
    if len(parts) < 2 * n_mat:
        raise ValueError(f"{name!r} lists fewer than {n_mat} names + fractions")
    names = parts[:n_mat]
    fracs = [float(x) for x in parts[n_mat:2 * n_mat]]
    return {
        "energy": energy,
        "n_materials": n_mat,
        "poly_ratio": int(poly_str) / 10.0,
        "thickness": f"oneL{geo_str}",
        "thickness_mm": int(geo_str),
        "materials": dict(zip(names, fracs)),
    }


# ---------------------------------------------------------------------------
# Namers
# ---------------------------------------------------------------------------
class CompactNamer:
    """The ``<beam>_<n>_<poly*10>_oneL<mm*100>_<names>_<fracs>`` scheme.

    Refuses (rather than silently mislabels) polymer fractions that are not a
    multiple of 0.1 and thicknesses that are not a multiple of 0.01 mm; use
    :class:`DescriptiveNamer` for those designs.
    """

    def name(self, spec):
        mix = spec.mixture
        poly_code = round(mix.polymer_fraction * 10)
        if abs(poly_code - mix.polymer_fraction * 10) > 1e-6:
            raise ValueError(
                f"CompactNamer encodes the polymer fraction in tenths; "
                f"{mix.polymer_fraction} would be mislabelled. "
                "Use Study(..., namer=attenupy.DescriptiveNamer())."
            )
        return (f"{spec.beam}_{mix.n_constituents}_{poly_code:02d}_"
                f"{mm_to_token(spec.thickness_mm)}_{mix.label}")

    @staticmethod
    def parse(name, polymer="C2H4"):
        """Rebuild a :class:`SimulationSpec` from a folder name."""
        info = parse_folder_name(name)
        return SimulationSpec(
            beam=info["energy"],
            mixture=Mixture(info["materials"], info["poly_ratio"], polymer),
            thickness_mm=info["thickness_mm"] / 100.0,
        )


class DescriptiveNamer:
    """Unambiguous names for arbitrary fractions/thicknesses, e.g.

    ``RQR6_080_n2_p25_t0.15mm_BaSO4-0.2500_Bi2O3-0.5000``
    """

    def __init__(self, decimals=4):
        self.decimals = decimals

    def name(self, spec):
        mix = spec.mixture
        parts = "_".join(f"{n}-{f:.{self.decimals}f}" for n, f in mix.fractions)
        poly = f"{mix.polymer_fraction * 100:g}"
        return (f"{spec.beam}_n{mix.n_constituents}_p{poly}_"
                f"t{spec.thickness_mm:g}mm_{parts}")
