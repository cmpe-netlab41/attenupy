"""Constituents: the named ingredients (compounds or elements) of a mixture.

A constituent pairs a *name* (what appears in folder names and CSV columns)
with a chemical *formula* (what determines the element weight fractions) and a
bulk *density*. Usually name and formula coincide (``"Bi2O3"``), but they don't
have to::

    import attenupy as atp
    atp.register_constituent("Gd2O3", density=7.41)
    atp.register_constituent("PU", density=1.2, formula="C3H8N2O")  # label != formula

Registries are plain objects, so a study can also use its own
:class:`ConstituentRegistry` instead of the process-wide default.
"""

from dataclasses import dataclass

from .chemistry import molar_mass, parse_formula
from .data import densities as _DEFAULT_DENSITIES


@dataclass(frozen=True)
class Constituent:
    name: str
    density: float          # g/cm**3
    formula: str = ""       # defaults to ``name``

    def __post_init__(self):
        if not self.formula:
            object.__setattr__(self, "formula", self.name)
        if not self.density or self.density <= 0:
            raise ValueError(f"Constituent {self.name!r}: density must be > 0")
        if "_" in self.name or " " in self.name:
            raise ValueError(
                f"Constituent name {self.name!r} may not contain '_' or spaces "
                "(it is used as a separator in folder names)"
            )
        parse_formula(self.formula)  # validate early

    @property
    def elements(self):
        return parse_formula(self.formula)

    @property
    def molar_mass(self):
        return molar_mass(self.formula)


class ConstituentRegistry:
    """A name -> :class:`Constituent` lookup with helpful errors."""

    def __init__(self, constituents=()):
        self._items = {}
        for c in constituents:
            self._items[c.name] = c

    @classmethod
    def from_densities(cls, densities):
        return cls(Constituent(name, d) for name, d in densities.items())

    def register(self, name, density, formula=None, *, overwrite=False):
        """Add a constituent. Returns it. Re-registering needs ``overwrite=True``."""
        c = Constituent(name, float(density), formula or "")
        existing = self._items.get(name)
        if existing is not None and existing != c and not overwrite:
            raise ValueError(
                f"Constituent {name!r} is already registered as {existing}; "
                "pass overwrite=True to replace it"
            )
        self._items[name] = c
        return c

    def resolve(self, item):
        """Return the :class:`Constituent` for a name (or pass one through)."""
        if isinstance(item, Constituent):
            return item
        try:
            return self._items[item]
        except KeyError:
            raise KeyError(
                f"Unknown constituent {item!r}. Register it first, e.g. "
                f"attenupy.register_constituent({item!r}, density=...)"
            ) from None

    __getitem__ = resolve

    def __contains__(self, name):
        return name in self._items

    def __iter__(self):
        return iter(self._items.values())

    def __len__(self):
        return len(self._items)

    def names(self):
        return sorted(self._items)

    def copy(self):
        return ConstituentRegistry(self._items.values())


DEFAULT_CONSTITUENTS = ConstituentRegistry.from_densities(_DEFAULT_DENSITIES)


def register_constituent(name, density, formula=None, *, overwrite=False):
    """Register a constituent in the default, process-wide registry."""
    return DEFAULT_CONSTITUENTS.register(name, density, formula, overwrite=overwrite)
