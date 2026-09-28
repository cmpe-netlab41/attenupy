"""Energy spectra for penEasy's BIGS source (the ``.spc`` file).

A spectrum source is any object with::

    write(path) -> None        # write a BIGS spectrum file at ``path``
    validate() -> None         # optional: raise if it cannot be produced

Built in:

* :class:`FileSpectrum` -- copy an existing ``.spc`` (user-provided or bundled).
* :class:`TabulatedSpectrum` -- write a spectrum from ``(energy, weight)`` arrays,
  e.g. the output of an X-ray tube model.
* :class:`SpekpySpectrum` -- compute an X-ray tube spectrum with SpekPy
  (optional dependency: ``pip install spekpy``).
* :func:`bundled_spectrum` -- the five beams used in the paper, shipped with the
  package (``N30_030``, ``RQR5_070``, ``RQR6_080``, ``RQR8_100``, ``x120kV``).

BIGS format, as penEasy reads it (``sourceBoxIsotropicGaussSpectrum.F``): one
``energy(eV) weight`` pair per line, the weight belonging to the bin from that
energy up to the next line's energy (sampled uniformly within it); weights are
relative (penEasy normalises them); lines with ``#`` in column 1 and blank lines
are skipped; the first line with a negative weight ends the list, and its energy
closes the last bin. :func:`check_bigs_spectrum` applies exactly these rules.

A spectrum known as *point* values f(E_k) (what most spectrum calculators
print) is not a histogram yet; :meth:`TabulatedSpectrum.from_points` converts it
the way the penEasy manual (Appendix B) prescribes.
"""

import os
import shutil
from dataclasses import dataclass, field

BUNDLED_DIR = os.path.join(os.path.dirname(__file__), "resources", "spectra")

# penEasy stops with "too many bins in spectrum" once the line count reaches
# NEMAX=32000 (sourceBoxIsotropicGaussSpectrum.F), so 31999 bins at most.
MAX_BINS = 31999


def read_spectrum(path):
    """Parse a BIGS ``.spc`` file into ``(energies_eV, weights)`` lists.

    The closing negative-probability line is not included.
    """
    energies, weights = [], []
    with open(path, "rb") as f:  # headers are often not valid UTF-8
        for raw in f:
            line = raw.decode("latin-1").strip()
            if not line or line.startswith("#"):
                continue
            e, w = (float(x) for x in line.split()[:2])
            if w < 0:
                break
            energies.append(e)
            weights.append(w)
    return energies, weights


def check_bigs_spectrum(path):
    """Check a ``.spc`` file the way penEasy's reader will; return the bin count.

    Raises :class:`ValueError` naming the offending line for anything that
    makes penEasy stop: a data line that isn't two numbers (e.g. a header
    without ``#`` in column 1), a negative or decreasing energy, more than
    :data:`MAX_BINS` bins, a missing negative-weight terminator, or no positive
    weight at all.
    """
    bins, total, last_e = 0, 0.0, None
    with open(path, "rb") as f:
        for lineno, raw in enumerate(f, 1):
            line = raw.decode("latin-1").rstrip("\r\n")[:80]  # penEasy reads a80
            if line.startswith("#") or not line.strip():
                continue
            try:
                e, w = (float(x) for x in line.replace(",", " ").split()[:2])
            except ValueError:
                raise ValueError(f"{path}:{lineno}: penEasy expects 'energy weight', got "
                                 f"{line.strip()!r} (comments need '#' in column 1)") from None
            if e < 0:
                raise ValueError(f"{path}:{lineno}: negative energy {e:g}")
            if last_e is not None and e < last_e:
                raise ValueError(f"{path}:{lineno}: decreasing energy {e:g} < {last_e:g}")
            last_e = e
            if w < 0:
                if bins == 0 or total <= 0:
                    raise ValueError(f"{path}: spectrum has no positive weight")
                return bins
            bins += 1
            total += w
            if bins > MAX_BINS:
                raise ValueError(f"{path}: more than {MAX_BINS} bins; penEasy's limit")
    raise ValueError(f"{path}: no closing line with a negative weight (e.g. "
                     f"'{last_e or 0:g}  -1'); penEasy would stop at end of file")


def write_bigs_spectrum(path, energies_eV, weights, *, header=None, end_energy_eV=None):
    """Write a BIGS spectrum file (a histogram).

    ``energies_eV`` are the lower bin edges, ``weights`` the (relative)
    probability of each bin, and ``end_energy_eV`` the upper edge of the last
    bin. By default the last bin is as wide as the one before it; a single
    energy with no ``end_energy_eV`` is a monoenergetic line.
    """
    if len(energies_eV) != len(weights) or not energies_eV:
        raise ValueError("energies and weights must be non-empty and equally long")
    if len(weights) > MAX_BINS:
        raise ValueError(f"{len(weights)} bins exceed penEasy's limit of {MAX_BINS}")
    if any(b <= a for a, b in zip(energies_eV, energies_eV[1:])):
        raise ValueError("energies must be strictly increasing")
    if any(w < 0 for w in weights):
        raise ValueError("weights must be >= 0")
    if end_energy_eV is not None:
        end = end_energy_eV
    elif len(energies_eV) > 1:
        end = 2 * energies_eV[-1] - energies_eV[-2]
    else:
        end = energies_eV[-1]
    if end < energies_eV[-1]:
        raise ValueError(f"end_energy_eV={end:g} is below the last bin edge")
    with open(path, "w") as f:
        for line in (header or "").splitlines():
            f.write((line if line.startswith("#") else f"# {line}") + "\n")
        f.write("# Energy(eV)  Probability(relative)\n")
        for e, w in zip(energies_eV, weights):
            f.write(f"{e:.8g}\t{w:.8g}\n")
        f.write(f"{end:.8g}\t-1\n")


@dataclass(frozen=True)
class FileSpectrum:
    """Copy an existing ``.spc`` file verbatim."""

    path: str

    def __post_init__(self):
        object.__setattr__(self, "path", os.path.abspath(self.path))

    def validate(self):
        if not os.path.isfile(self.path):
            raise FileNotFoundError(f"Spectrum file not found: {self.path}")
        check_bigs_spectrum(self.path)

    def write(self, path):
        shutil.copyfile(self.path, path)

    def read(self):
        return read_spectrum(self.path)


@dataclass(frozen=True)
class TabulatedSpectrum:
    """A histogram spectrum: lower bin edges (eV) and relative bin weights.

    For point values of a continuous spectrum use :meth:`from_points`.
    """

    energies_eV: tuple
    weights: tuple
    header: str = None
    end_energy_eV: float = None

    def __post_init__(self):
        object.__setattr__(self, "energies_eV", tuple(float(e) for e in self.energies_eV))
        object.__setattr__(self, "weights", tuple(float(w) for w in self.weights))

    @classmethod
    def from_points(cls, energies_eV, values, **kwargs):
        """Histogram from point values f(E_k) of a continuous spectrum.

        Each bin [E_k, E_k+1] gets the trapezoid ``(f_k + f_k+1)/2 * (E_k+1 - E_k)``
        (penEasy manual, Appendix B), and the last point closes the last bin.
        """
        e, f = [float(x) for x in energies_eV], [float(x) for x in values]
        if len(e) != len(f) or len(e) < 2:
            raise ValueError("need at least two points, as many energies as values")
        weights = [(f[k] + f[k + 1]) / 2 * (e[k + 1] - e[k]) for k in range(len(e) - 1)]
        return cls(e[:-1], weights, end_energy_eV=e[-1], **kwargs)

    def validate(self):
        if len(self.energies_eV) != len(self.weights) or not self.weights:
            raise ValueError("energies and weights must be non-empty and equally long")

    def write(self, path):
        write_bigs_spectrum(path, self.energies_eV, self.weights,
                            header=self.header, end_energy_eV=self.end_energy_eV)

    def __repr__(self):
        e = self.energies_eV
        return f"TabulatedSpectrum({len(e)} bins, {e[0]:g}-{e[-1]:g} eV)"


@dataclass(frozen=True)
class SpekpySpectrum:
    """An X-ray tube spectrum computed with SpekPy (``pip install spekpy``).

        SpekpySpectrum(kvp=70, th=12, filters=[("Al", 2.8), ("Be", 3.0)])

    ``th`` is the anode angle (deg), ``dk`` the bin width (keV), ``z`` the
    distance (cm); ``filters`` are ``(material, mm)`` pairs applied in order;
    ``spek_kwargs`` go to ``spekpy.Spek`` unchanged (e.g. ``targ="Mo"``).

    Each output row is a SpekPy bin: its lower edge in eV and its fluence
    (photons/cm2/mAs/keV, unrounded; penEasy only uses the relative values),
    closed by ``<kVp in eV>  -1``. The header records the parameters, the HVLs
    and the mean energy.
    """

    kvp: float
    th: float = 12.0
    filters: tuple = ()
    dk: float = 1.0
    z: float = 100.0
    spek_kwargs: dict = field(default_factory=dict, hash=False)

    def __post_init__(self):
        object.__setattr__(self, "kvp", float(self.kvp))
        object.__setattr__(self, "filters",
                           tuple((str(m), float(mm)) for m, mm in self.filters))

    def spek(self):
        """The configured ``spekpy.Spek`` object."""
        try:
            import spekpy
        except ImportError:
            raise ImportError("SpekpySpectrum needs SpekPy: pip install spekpy") from None
        s = spekpy.Spek(kvp=self.kvp, th=self.th, dk=self.dk, z=self.z, **self.spek_kwargs)
        for material, mm in self.filters:
            s.filter(material, mm)
        return s

    def validate(self):
        self.spek()

    def table(self):
        """``(lower_edges_eV, fluence, upper_edge_eV)`` of the spectrum."""
        k, f = self.spek().get_spectrum()          # bin centres (keV), fluence
        edges = [round((kc - self.dk / 2) * 1000, 6) for kc in k]
        return edges, [float(x) for x in f], round((k[-1] + self.dk / 2) * 1000, 6)

    def hvl(self):
        """``(HVL1, HVL2)`` in mm Al, e.g. to check against IEC 61267."""
        s = self.spek()
        return s.get_hvl1(), s.get_hvl2()

    def write(self, path):
        s = self.spek()
        edges, fluence, end = self.table()
        mean = sum((e / 1000 + self.dk / 2) * w for e, w in zip(edges, fluence)) / sum(fluence)
        filters = ", ".join(f"{m} {mm:g} mm" for m, mm in self.filters) or "none"
        header = (f"X-ray spectrum from SpekPy: {self.kvp:g} kVp, anode angle {self.th:g} deg, "
                  f"distance {self.z:g} cm, bin width {self.dk:g} keV\n"
                  f"Filtration: {filters}\n"
                  f"HVL1 = {s.get_hvl1():.4f} mm Al, HVL2 = {s.get_hvl2():.4f} mm Al, "
                  f"mean energy = {mean:.3f} keV\n"
                  "Rows: bin lower edge (eV), fluence (photons/cm2/mAs/keV)")
        if self.spek_kwargs:
            header += f"\nSpek options: {self.spek_kwargs}"
        write_bigs_spectrum(path, edges, fluence, header=header, end_energy_eV=end)


def bundled_spectra():
    """Names of the spectra shipped with attenupy."""
    return sorted(f[:-4] for f in os.listdir(BUNDLED_DIR) if f.endswith(".spc"))


def bundled_spectrum(name):
    """A :class:`FileSpectrum` for one of the shipped spectra, e.g. ``"RQR6_080"``."""
    path = os.path.join(BUNDLED_DIR, f"{name}.spc")
    if not os.path.isfile(path):
        raise KeyError(f"No bundled spectrum {name!r}; available: {bundled_spectra()}")
    return FileSpectrum(path)
