"""Reference data: elements, default constituent densities, beam tables.

Atomic weights are IUPAC conventional values for Z = 1..99 (the range
PENELOPE supports), with a small override table for a few elements where a
slightly different value is kept so that element weight fractions -- and
therefore every ``mat.mat`` -- stay consistent. Densities are in g/cm**3,
atomic weights in g/mol.
"""

# ---------------------------------------------------------------------------
# Elements
# ---------------------------------------------------------------------------
# (symbol, Z, standard atomic weight)
_ELEMENTS = [
    ("H", 1, 1.008), ("He", 2, 4.002602), ("Li", 3, 6.94), ("Be", 4, 9.0121831),
    ("B", 5, 10.81), ("C", 6, 12.011), ("N", 7, 14.007), ("O", 8, 15.999),
    ("F", 9, 18.998403163), ("Ne", 10, 20.1797), ("Na", 11, 22.98976928),
    ("Mg", 12, 24.305), ("Al", 13, 26.9815385), ("Si", 14, 28.085),
    ("P", 15, 30.973761998), ("S", 16, 32.06), ("Cl", 17, 35.45),
    ("Ar", 18, 39.948), ("K", 19, 39.0983), ("Ca", 20, 40.078),
    ("Sc", 21, 44.955908), ("Ti", 22, 47.867), ("V", 23, 50.9415),
    ("Cr", 24, 51.9961), ("Mn", 25, 54.938044), ("Fe", 26, 55.845),
    ("Co", 27, 58.933194), ("Ni", 28, 58.6934), ("Cu", 29, 63.546),
    ("Zn", 30, 65.38), ("Ga", 31, 69.723), ("Ge", 32, 72.630),
    ("As", 33, 74.921595), ("Se", 34, 78.971), ("Br", 35, 79.904),
    ("Kr", 36, 83.798), ("Rb", 37, 85.4678), ("Sr", 38, 87.62),
    ("Y", 39, 88.90584), ("Zr", 40, 91.224), ("Nb", 41, 92.90637),
    ("Mo", 42, 95.95), ("Tc", 43, 98.0), ("Ru", 44, 101.07),
    ("Rh", 45, 102.90550), ("Pd", 46, 106.42), ("Ag", 47, 107.8682),
    ("Cd", 48, 112.414), ("In", 49, 114.818), ("Sn", 50, 118.710),
    ("Sb", 51, 121.760), ("Te", 52, 127.60), ("I", 53, 126.90447),
    ("Xe", 54, 131.293), ("Cs", 55, 132.90545196), ("Ba", 56, 137.327),
    ("La", 57, 138.90547), ("Ce", 58, 140.116), ("Pr", 59, 140.90766),
    ("Nd", 60, 144.242), ("Pm", 61, 145.0), ("Sm", 62, 150.36),
    ("Eu", 63, 151.964), ("Gd", 64, 157.25), ("Tb", 65, 158.92535),
    ("Dy", 66, 162.500), ("Ho", 67, 164.93033), ("Er", 68, 167.259),
    ("Tm", 69, 168.93422), ("Yb", 70, 173.045), ("Lu", 71, 174.9668),
    ("Hf", 72, 178.49), ("Ta", 73, 180.94788), ("W", 74, 183.84),
    ("Re", 75, 186.207), ("Os", 76, 190.23), ("Ir", 77, 192.217),
    ("Pt", 78, 195.084), ("Au", 79, 196.966569), ("Hg", 80, 200.592),
    ("Tl", 81, 204.38), ("Pb", 82, 207.2), ("Bi", 83, 208.98040),
    ("Po", 84, 209.0), ("At", 85, 210.0), ("Rn", 86, 222.0),
    ("Fr", 87, 223.0), ("Ra", 88, 226.0), ("Ac", 89, 227.0),
    ("Th", 90, 232.0377), ("Pa", 91, 231.03588), ("U", 92, 238.02891),
    ("Np", 93, 237.0), ("Pu", 94, 244.0), ("Am", 95, 243.0),
    ("Cm", 96, 247.0), ("Bk", 97, 247.0), ("Cf", 98, 251.0),
    ("Es", 99, 252.0),
]

# Override values for a few elements (see the module docstring).
_ATOMIC_WEIGHT_OVERRIDES = {
    "Bi": 208.9804, "W": 183.84, "Sb": 121.76, "Sn": 118.71,
    "Cu": 63.546, "B": 10.81, "Ba": 137.327,
    "O": 15.999, "S": 32.065, "C": 12.0107, "H": 1.00794, "N": 14.0067,
    "Pb": 207.231, "Gd": 157.25,
}

atomic_numbers = {sym: z for sym, z, _ in _ELEMENTS}
atomic_weights = {sym: w for sym, _, w in _ELEMENTS}
atomic_weights.update(_ATOMIC_WEIGHT_OVERRIDES)

# ---------------------------------------------------------------------------
# Constituents (compounds / elements that can go into a mixture)
# ---------------------------------------------------------------------------
# Polymer (matrix) phase. PENELOPE/penEasy treats it as the compound C2H4.
POLYMER_FORMULA = "C2H4"
POLYMER_DENSITY = 0.95  # g/cm**3

# Densities of the shielding constituents shipped by default. Anything else is
# added at run time with ``attenupy.register_constituent(name, density=...)``.
densities = {
    "Bi2O3": 8.9018, "BaSO4": 4.48, "BN": 2.1,
    "WO3": 7.1596, "B4C": 2.52,
    "H": 0.00008988, "C": 2.266, "N": 0.001251,
    "O": 0.001429, "S": 1.96,
    "Bi": 9.747, "W": 19.254, "Sb": 6.697,
    "Sn": 7.265, "Cu": 8.96, "B": 2.37,
    "Ba": 3.59, "Pb": 11.34, "Gd": 7.90,
    POLYMER_FORMULA: POLYMER_DENSITY,
}

# PENELOPE pdcompos.pen material ID for dry air (used for Aire.mat).
AIR_MATERIAL_ID = 104

# ---------------------------------------------------------------------------
# Built-in beams (registered into attenupy.beams.DEFAULT_BEAMS)
# ---------------------------------------------------------------------------
ENERGY_TO_FOLDER = {
    "N30_030": "30kV",
    "RQR5_070": "70kV",
    "RQR6_080": "80kV",
    "RQR8_100": "100kV",
    "x120kV": "120kV",
}

# Energy deposited in the air detector with the slab replaced by air (eV per
# history), i.e. I0 for the default slab geometry. Keyed by subfolder name.
AIRE_ENERGY = {
    "30kV": 1.92544e00,
    "70kV": 9.35590e-01,
    "80kV": 8.71549e-01,
    "100kV": 7.73834e-01,
    "120kV": 6.30952e-01,
}
