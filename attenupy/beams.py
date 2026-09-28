"""Beams: an X-ray tube setting = a key, a tube voltage and a spectrum.

The five beams of the paper are pre-registered
(with their spectra bundled), so they work out of the box. New ones::

    import attenupy as atp
    atp.register_beam("RQR9_120", kvp=120, spectrum="spectra/RQR9_120.spc")
    atp.register_beam("W60", kvp=60, spectrum=atp.TabulatedSpectrum(E, w))

``reference_energy`` is I0: the detector signal (eV/history) with the slab
replaced by air, used to turn a tally into a transmission ratio. It is only
valid for the geometry it was measured with; leave it ``None`` for new beams
and let the study measure it with an air reference run
(:meth:`attenupy.Study.prepare_references`).
"""

import re
from dataclasses import dataclass

from .data import AIRE_ENERGY, ENERGY_TO_FOLDER
from .spectra import FileSpectrum, bundled_spectrum

_KEY_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9.+-]*(?:_[A-Za-z0-9.+-]+)*$")


@dataclass(frozen=True)
class Beam:
    key: str                        # "RQR6_080": folder-name prefix and .spc stem
    kvp: float                      # tube voltage (kV)
    spectrum: object                # a spectrum source (see attenupy.spectra)
    folder: str = ""                # sub-directory of the study; default "<kvp>kV"
    reference_energy: float = None  # I0 in eV/history, or None to measure it

    def __post_init__(self):
        if not _KEY_RE.match(self.key) or len(self.key) > 60:
            raise ValueError(f"Invalid beam key {self.key!r} (letters, digits, '_', '.', '-')")
        if not self.folder:
            object.__setattr__(self, "folder", f"{self.kvp:g}kV")

    @property
    def spectrum_filename(self):
        return f"{self.key}.spc"


class BeamRegistry:
    """Key/alias -> :class:`Beam`. Aliases: the kV number and folder name."""

    def __init__(self, beams=()):
        self._beams = {}
        self._aliases = {}
        for b in beams:
            self.add(b)

    def add(self, beam, *, aliases=(), overwrite=False):
        if beam.key in self._beams and not overwrite and self._beams[beam.key] != beam:
            raise ValueError(f"Beam {beam.key!r} already registered; pass overwrite=True")
        self._beams[beam.key] = beam
        for alias in (f"{beam.kvp:g}", beam.folder, *aliases):
            if alias == beam.key:
                continue
            owner = self._aliases.get(alias, beam.key)
            # An alias shared by two beams becomes ambiguous (None) rather than
            # silently pointing at whichever was registered last.
            self._aliases[alias] = beam.key if owner == beam.key else None
        return beam

    def register(self, key, kvp, spectrum=None, *, folder="", reference_energy=None,
                 aliases=(), overwrite=False):
        """Create and add a :class:`Beam`.

        ``spectrum`` may be a spectrum source, a path to a ``.spc`` file, or
        ``None`` to use the bundled spectrum called ``key``.
        """
        if spectrum is None:
            spectrum = bundled_spectrum(key)
        elif isinstance(spectrum, str):
            spectrum = FileSpectrum(spectrum)
        beam = Beam(key, float(kvp), spectrum, folder, reference_energy)
        return self.add(beam, aliases=aliases, overwrite=overwrite)

    def resolve(self, item):
        """Beam for a :class:`Beam`, key (``"RQR6_080"``) or alias (``"80"``, ``"80kV"``)."""
        if isinstance(item, Beam):
            return item
        item = str(item)
        if item in self._beams:
            return self._beams[item]
        if item in self._aliases:
            key = self._aliases[item]
            if key is None:
                matches = [b.key for b in self._beams.values()
                           if item in (f"{b.kvp:g}", b.folder)]
                raise KeyError(f"Beam alias {item!r} is ambiguous ({matches}); use the key")
            return self._beams[key]
        raise KeyError(
            f"Unknown beam {item!r}. Known: {sorted(self._beams)}. Register one with "
            f"attenupy.register_beam({item!r}, kvp=..., spectrum=...)"
        )

    __getitem__ = resolve

    def by_folder(self, folder):
        for b in self._beams.values():
            if b.folder == folder:
                return b
        raise KeyError(f"No beam uses folder {folder!r}")

    def __contains__(self, item):
        try:
            self.resolve(item)
            return True
        except KeyError:
            return False

    def __iter__(self):
        return iter(self._beams.values())

    def __len__(self):
        return len(self._beams)

    def copy(self):
        reg = BeamRegistry()
        reg._beams = dict(self._beams)
        reg._aliases = dict(self._aliases)
        return reg


_PAPER_KVP = {"N30_030": 30, "RQR5_070": 70, "RQR6_080": 80, "RQR8_100": 100, "x120kV": 120}

DEFAULT_BEAMS = BeamRegistry()
for _key, _kvp in _PAPER_KVP.items():
    DEFAULT_BEAMS.register(_key, _kvp, folder=ENERGY_TO_FOLDER[_key],
                           reference_energy=AIRE_ENERGY[ENERGY_TO_FOLDER[_key]])


def register_beam(key, kvp, spectrum=None, **kwargs):
    """Register a beam in the default, process-wide registry."""
    return DEFAULT_BEAMS.register(key, kvp, spectrum, **kwargs)
