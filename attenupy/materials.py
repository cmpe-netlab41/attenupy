"""Material files (``mat.mat``, ``Aire.mat``) via PENELOPE's ``material.exe``.

A material builder is any object with::

    build(mixture, dest_path, registry) -> None   # write the slab material
    build_air(dest_path) -> None                  # write dry air

:class:`PenelopeMaterials` implements it on top of
:class:`PenelopeMaterialTool`, a small driver for the interactive
``material.exe``. The driver runs the tool in a private temporary directory
(with a ``pdfiles`` link to the database), so it

* never hits material.exe's 80-character output-path limit,
* works when the PENELOPE install is read-only,
* can run several instances concurrently.
"""

import os
import shutil
import subprocess
import tempfile

from .constituents import DEFAULT_CONSTITUENTS
from .data import AIR_MATERIAL_ID, atomic_numbers

# Answers to material.exe's interactive prompts.
_FROM_KEYBOARD = "1"
_FROM_DATABASE = "2"
_BY_WEIGHT = "2"
_YES, _NO = "1", "2"
_MAX_NAME = 62  # material.exe reads the name as A62


class MaterialToolError(RuntimeError):
    pass


def composition_input(name, z_fractions, density, output_name, mean_excitation_eV=None):
    """stdin for material.exe describing a material by element weight fractions.

    ``z_fractions`` maps atomic number -> weight fraction (any normalisation).
    """
    lines = [_FROM_KEYBOARD, name[:_MAX_NAME], str(len(z_fractions))]
    if len(z_fractions) > 1:
        lines.append(_BY_WEIGHT)
    for z in sorted(z_fractions):
        lines.append(f"{z} {z_fractions[z]:.6f}")
    if mean_excitation_eV is None:
        lines.append(_NO)
    else:
        lines += [_YES, f"{mean_excitation_eV:.6f}"]
    lines += [f"{density:.6f}", _NO, output_name]
    return "\n".join(lines) + "\n"


def database_input(material_id, output_name):
    """stdin for material.exe loading entry ``material_id`` of pdcompos.pen."""
    return "\n".join([_FROM_DATABASE, str(int(material_id)), _NO, output_name]) + "\n"


def mixture_input(mixture, output_name, registry=None):
    """``(stdin, density)`` for a :class:`~attenupy.spec.Mixture`.

    Element fractions are listed in element-symbol order.
    """
    fractions = mixture.element_fractions(registry)
    density = mixture.density(registry)
    lines = [_FROM_KEYBOARD, mixture.label[:_MAX_NAME], str(len(fractions))]
    if len(fractions) > 1:
        lines.append(_BY_WEIGHT)
    for el in sorted(fractions):
        lines.append(f"{atomic_numbers[el]} {fractions[el]:.6f}")
    lines += [_NO, f"{density:.6f}", _NO, output_name]
    return "\n".join(lines), density


class PenelopeMaterialTool:
    """Drive ``material.exe`` non-interactively."""

    def __init__(self, exe, *, timeout=900):
        self.exe = os.path.abspath(os.path.expanduser(exe))
        self.pendbase = os.path.dirname(self.exe)
        self.pdfiles = os.path.join(self.pendbase, "pdfiles")
        self.timeout = timeout

    def run(self, stdin, dest, output_name="out.mat"):
        """Feed ``stdin`` (which must name ``output_name`` last) and move the
        resulting file to ``dest``. Returns material.exe's stdout."""
        workdir = tempfile.mkdtemp(prefix="attenupy-mat-")
        try:
            try:
                os.symlink(self.pdfiles, os.path.join(workdir, "pdfiles"))
                cwd, out = workdir, output_name
            except OSError:  # no symlinks: fall back to the pendbase dir itself
                cwd = self.pendbase
                out = f"ps{os.getpid()}_{os.path.basename(workdir)[-6:]}.mat"
                stdin = stdin.replace(output_name, out)
            proc = subprocess.run(
                [self.exe], input=stdin, cwd=cwd, capture_output=True, text=True,
                timeout=self.timeout,
            )
            produced = os.path.join(cwd, out)
            if proc.returncode != 0 or not os.path.isfile(produced) \
                    or os.path.getsize(produced) == 0:
                tail = "\n".join((proc.stdout + proc.stderr).splitlines()[-15:])
                raise MaterialToolError(
                    f"material.exe failed (exit {proc.returncode}) for {dest}:\n{tail}"
                )
            os.makedirs(os.path.dirname(os.path.abspath(dest)), exist_ok=True)
            shutil.move(produced, dest)
            return proc.stdout
        finally:
            shutil.rmtree(workdir, ignore_errors=True)

    def from_composition(self, name, z_fractions, density, dest, mean_excitation_eV=None):
        return self.run(composition_input(name, z_fractions, density, "out.mat",
                                          mean_excitation_eV), dest)

    def from_database(self, material_id, dest):
        return self.run(database_input(material_id, "out.mat"), dest)


class PenelopeMaterials:
    """Material builder backed by ``material.exe``."""

    def __init__(self, tool):
        if not isinstance(tool, PenelopeMaterialTool):
            tool = PenelopeMaterialTool(tool)
        self.tool = tool

    def build(self, mixture, dest, registry=None):
        stdin, _ = mixture_input(mixture, "out.mat", registry or DEFAULT_CONSTITUENTS)
        self.tool.run(stdin + "\n", dest)

    def build_air(self, dest):
        self.tool.from_database(AIR_MATERIAL_ID, dest)
