"""Describe a study design in a JSON or TOML file instead of Python.

Used by ``attenupy prepare ROOT --design design.toml``. Example::

    [constituents]            # extra constituents: name = density or a table
    Gd2O3 = 7.41
    PU = { density = 1.2, formula = "C3H8N2O" }

    [[beams]]                 # extra beams (the paper's five are built in)
    key = "W60"
    kvp = 60
    spectrum = "spectra/W60.spc"          # a .spc file, or computed with SpekPy:

    [[beams]]
    key = "RQR5_spekpy"
    kvp = 70
    spectrum = { type = "spekpy", th = 12, filters = [["Al", 2.8], ["Be", 3.0]] }

    [[design]]                # the paper's training design
    type = "random"
    n = 300
    beams = ["RQR6_080", "RQR8_100", "x120kV"]
    constituents = ["Bi2O3", "BaSO4", "BN", "B4C", "WO3", "Sb", "Sn", "Bi", "W", "Cu", "B", "Ba"]
    n_constituents = [1, 2, 3]
    polymer_fractions = [0.2, 0.3, 0.4]
    thicknesses_mm = [0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.05, 1.2, 1.35, 1.5]
    seed = 0
    splitter = { type = "dirichlet", alpha = 1.0 }

    [[design]]                # a held-out thickness curve
    type = "sweep"
    beams = ["RQR6_080"]
    thicknesses_mm = { start = 0.15, stop = 1.5, step = 0.15 }
    mixtures = [ { fractions = { Bi2O3 = 0.4, BaSO4 = 0.2 }, polymer_fraction = 0.4 } ]

Custom samplers become available to design files with :func:`register_sampler`.
"""

import json
import os

from .beams import DEFAULT_BEAMS
from .constituents import DEFAULT_CONSTITUENTS
from .sampling import (Chain, CSVSampler, Dirichlet, EqualSplit, GridSampler,
                       RandomMixtures, RandomSampler, ThicknessSweep, thickness_range)
from .spec import Mixture
from .spectra import FileSpectrum, SpekpySpectrum

SAMPLERS = {}
SPLITTERS = {"dirichlet": Dirichlet, "equal": EqualSplit}
SPECTRA = {"spekpy": SpekpySpectrum, "file": FileSpectrum}


def register_sampler(name):
    """Decorator: make a sampler class usable as ``type = "<name>"``.

    The class is constructed with the remaining keys of its ``[[design]]``
    table as keyword arguments.
    """
    def deco(cls):
        SAMPLERS[name] = cls
        return cls
    return deco


def _thicknesses(value):
    if isinstance(value, dict):
        return thickness_range(value["start"], value["stop"], value["step"])
    return list(value)


def _splitter(value):
    if value is None or not isinstance(value, dict):
        return value
    value = dict(value)
    return SPLITTERS[value.pop("type", "dirichlet")](**value)


def _mixtures(value):
    if isinstance(value, dict) and "random" in value:  # {random = {...}, n = 6, seed = 1}
        spec = dict(value["random"])
        spec["splitter"] = _splitter(spec.get("splitter"))
        return RandomMixtures(**spec).sample(value["n"], seed=value.get("seed"))
    return [Mixture.from_dict(m) for m in value]


def _random(**kw):
    kw["thicknesses_mm"] = _thicknesses(kw["thicknesses_mm"])
    if "splitter" in kw:
        kw["splitter"] = _splitter(kw["splitter"])
    return RandomSampler(**kw)


def _grid(**kw):
    return GridSampler(kw.pop("beams"), _mixtures(kw.pop("mixtures")),
                       _thicknesses(kw.pop("thicknesses_mm")), **kw)


def _sweep(**kw):
    return ThicknessSweep(_mixtures(kw.pop("mixtures")),
                          _thicknesses(kw.pop("thicknesses_mm")), kw.pop("beams"), **kw)


SAMPLERS.update({"random": _random, "grid": _grid, "sweep": _sweep, "csv": CSVSampler})


def load_design_file(path):
    """Read a ``.json`` or ``.toml`` design file into a dict."""
    if path.endswith(".toml"):
        try:
            import tomllib
        except ImportError:  # Python < 3.11
            try:
                import tomli as tomllib
            except ImportError:
                raise RuntimeError("Reading TOML needs Python >= 3.11 or `pip install tomli`; "
                                   "or write the design as JSON") from None
        with open(path, "rb") as f:
            return tomllib.load(f)
    with open(path) as f:
        return json.load(f)


def design_from_dict(config, *, base_dir=".", constituents=None, beams=None):
    """Register the file's constituents/beams and return one combined sampler."""
    constituents = constituents or DEFAULT_CONSTITUENTS
    beams = beams or DEFAULT_BEAMS
    for name, value in config.get("constituents", {}).items():
        if isinstance(value, dict):
            constituents.register(name, value["density"], value.get("formula"), overwrite=True)
        else:
            constituents.register(name, value, overwrite=True)
    for b in config.get("beams", []):
        b = dict(b)
        spectrum = b.pop("spectrum", None)
        if isinstance(spectrum, dict):  # { type = "spekpy", ... }; kvp defaults to the beam's
            spectrum = dict(spectrum)
            kind = spectrum.pop("type")
            if kind not in SPECTRA:
                raise ValueError(f"beam {b.get('key')!r}: unknown spectrum type {kind!r}; "
                                 f"known: {sorted(SPECTRA)}")
            if kind == "spekpy":
                spectrum.setdefault("kvp", b["kvp"])
            elif not os.path.isabs(spectrum.get("path", "/")):
                spectrum["path"] = os.path.join(base_dir, spectrum["path"])
            spectrum = SPECTRA[kind](**spectrum)
        elif isinstance(spectrum, str) and not os.path.isabs(spectrum):
            spectrum = os.path.join(base_dir, spectrum)
        beams.register(b.pop("key"), b.pop("kvp"), spectrum, overwrite=True, **b)

    parts = []
    for i, entry in enumerate(config.get("design", [])):
        entry = dict(entry)
        kind = entry.pop("type", None)
        if kind not in SAMPLERS:
            raise ValueError(f"design #{i}: unknown type {kind!r}; known: {sorted(SAMPLERS)}")
        if kind == "csv" and not os.path.isabs(entry.get("path", "/")):
            entry["path"] = os.path.join(base_dir, entry["path"])
        parts.append(SAMPLERS[kind](**entry))
    if not parts:
        raise ValueError("The design file has no [[design]] entries")
    return parts[0] if len(parts) == 1 else Chain(*parts)


def load_design(path, **kwargs):
    """Sampler described by the design file at ``path``."""
    return design_from_dict(load_design_file(path),
                            base_dir=os.path.dirname(os.path.abspath(path)), **kwargs)
