"""Chemistry helpers: formula parsing, element weight fractions, mixture density.

The formula parser raises on unknown tokens rather than silently skipping
them, and understands parentheses, e.g. ``"Mg(OH)2"`` or ``"Ca10(PO4)6(OH)2"``.
"""

import re

from .data import (
    POLYMER_DENSITY,
    POLYMER_FORMULA,
    atomic_weights,
    densities as DEFAULT_DENSITIES,
)

_TOKEN = re.compile(r"([A-Z][a-z]?|\(|\)|\d+(?:\.\d+)?)")


def parse_formula(formula):
    """Break a chemical formula into ``{element: count}``.

    ``"Bi2O3"`` -> ``{"Bi": 2, "O": 3}``; ``"Mg(OH)2"`` -> ``{"Mg": 1, "O": 2, "H": 2}``.
    Fractional counts (``"Fe0.5Ni0.5"``) are allowed. Raises :class:`ValueError`
    for anything that is not a formula made of known element symbols.
    """
    tokens = _TOKEN.findall(formula)
    if "".join(tokens) != formula or not formula:
        raise ValueError(f"Cannot parse chemical formula {formula!r}")

    stack = [{}]
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        i += 1
        if tok == "(":
            stack.append({})
            continue
        if tok[0].isdigit():
            raise ValueError(f"Unexpected number in formula {formula!r}")
        mult = 1
        if i < len(tokens) and tokens[i][0].isdigit():
            mult = float(tokens[i])
            mult = int(mult) if mult.is_integer() else mult
            i += 1
        if tok == ")":
            if len(stack) == 1:
                raise ValueError(f"Unbalanced ')' in formula {formula!r}")
            for el, c in stack.pop().items():
                stack[-1][el] = stack[-1].get(el, 0) + c * mult
        else:
            if tok not in atomic_weights:
                raise ValueError(f"Unknown element {tok!r} in formula {formula!r}")
            stack[-1][tok] = stack[-1].get(tok, 0) + mult
    if len(stack) != 1:
        raise ValueError(f"Unbalanced '(' in formula {formula!r}")
    return stack[0]


def molar_mass(formula):
    """Molar mass (g/mol) of ``formula``."""
    return sum(atomic_weights[el] * c for el, c in parse_formula(formula).items())


def element_weight_fractions(mixture):
    """Element-level weight fractions for a mixture of compounds.

    ``mixture`` is an iterable of ``(formula, weight)`` pairs, where ``weight``
    is the weight fraction of that compound in the mixture; PENELOPE normalises
    internally, so the values are used as-is.

    Returns ``{element: summed_weight_fraction}``.
    """
    overall = {}
    for formula, weight in mixture:
        counts = parse_formula(formula)
        mol_weight = sum(atomic_weights[el] * c for el, c in counts.items())
        for el, c in counts.items():
            frac_e = (atomic_weights[el] * c) / mol_weight
            overall[el] = overall.get(el, 0.0) + weight * frac_e
    return overall


def mixture_density(mixture, densities_dict=None):
    """Inverse-rule-of-mixtures density (g/cm**3) for ``(formula, weight)`` pairs.

    The polymer phase always uses :data:`POLYMER_DENSITY` regardless of what the
    densities table says, matching the original behaviour.
    """
    if densities_dict is None:
        densities_dict = DEFAULT_DENSITIES
    inv = 0.0
    for formula, weight in mixture:
        d = POLYMER_DENSITY if formula == POLYMER_FORMULA else densities_dict.get(formula)
        if d is None or d <= 0:
            raise ValueError(f"Missing or invalid density for {formula!r}")
        inv += weight / d
    return 1.0 / inv


def build_mixture(materials, poly_ratio):
    """Assemble the full ``(formula, weight)`` list, polymer included.

    ``materials`` maps non-polymer compound -> weight fraction. The polymer is
    prepended (as in the original) only when ``poly_ratio > 0``.
    """
    mixture = []
    if poly_ratio > 0.0:
        mixture.append((POLYMER_FORMULA, float(poly_ratio)))
    for formula, weight in materials.items():
        mixture.append((formula, float(weight)))
    return mixture
