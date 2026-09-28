"""The two user-supplied PENELOPE executables, and checks that they will run.

attenupy compiles nothing. You bring, built for the machine that will run
them:

* ``penEasy.x``    -- the penEasy main program;
* ``material.exe`` -- PENELOPE's material generator, sitting in its
  ``pendbase`` directory next to the ``pdfiles/`` database.

Point attenupy at them explicitly, or through environment variables::

    tc = Toolchain(peneasy="~/pen/penEasy.x", material="~/pen/pendbase/material.exe")
    tc = Toolchain.from_env()   # ATTENUPY_PENEASY / ATTENUPY_MATERIAL

``attenupy doctor`` (or ``Toolchain.problems()``) reports what is missing,
including executables built for a different OS/CPU than the current one.
"""

import os
import platform
import struct
from dataclasses import dataclass

ENV_PENEASY = "ATTENUPY_PENEASY"
ENV_MATERIAL = "ATTENUPY_MATERIAL"

_ELF_MACHINES = {0x3E: "x86_64", 0xB7: "arm64", 0x03: "x86"}
_MACHO_CPUS = {0x01000007: "x86_64", 0x0100000C: "arm64"}


def binary_platform(path):
    """Best-effort ``(os, arch)`` of an executable, e.g. ``("Linux", "x86_64")``.

    Returns ``None`` if the format is not recognised (e.g. a shell script).
    """
    try:
        with open(path, "rb") as f:
            head = f.read(20)
    except OSError:
        return None
    if head[:4] == b"\x7fELF":
        endian = "<" if head[5] == 1 else ">"
        (machine,) = struct.unpack(endian + "H", head[18:20])
        return ("Linux", _ELF_MACHINES.get(machine, hex(machine)))
    if head[:4] in (b"\xcf\xfa\xed\xfe", b"\xce\xfa\xed\xfe"):
        (cpu,) = struct.unpack("<I", head[4:8])
        return ("Darwin", _MACHO_CPUS.get(cpu, hex(cpu)))
    if head[:2] == b"MZ":
        return ("Windows", "?")
    return None


def _host_platform():
    machine = platform.machine().lower()
    arch = {"amd64": "x86_64", "aarch64": "arm64"}.get(machine, machine)
    return (platform.system(), arch)


def _check_executable(label, path):
    if not path:
        return [f"{label}: not configured"]
    if not os.path.isfile(path):
        return [f"{label}: {path} does not exist"]
    problems = []
    if not os.access(path, os.X_OK):
        problems.append(f"{label}: {path} is not executable (chmod +x)")
    found, host = binary_platform(path), _host_platform()
    if found and found[0] != "Windows" and found != host:
        problems.append(
            f"{label}: {path} is built for {found[0]}/{found[1]} but this machine is "
            f"{host[0]}/{host[1]}; it will only run on a matching machine"
        )
    return problems


def _expand(path):
    return os.path.abspath(os.path.expanduser(path)) if path else None


@dataclass
class Toolchain:
    peneasy: str = None    # path to penEasy.x
    material: str = None   # path to material.exe (inside pendbase/)

    def __post_init__(self):
        self.peneasy = _expand(self.peneasy)
        self.material = _expand(self.material)

    @classmethod
    def from_env(cls, peneasy=None, material=None):
        """Explicit arguments win; otherwise read the environment variables."""
        return cls(peneasy or os.environ.get(ENV_PENEASY),
                   material or os.environ.get(ENV_MATERIAL))

    @property
    def pdfiles(self):
        return os.path.join(os.path.dirname(self.material), "pdfiles") if self.material else None

    def peneasy_problems(self):
        return _check_executable("penEasy.x", self.peneasy)

    def material_problems(self):
        problems = _check_executable("material.exe", self.material)
        if self.material and not os.path.isfile(os.path.join(self.pdfiles, "pdcompos.pen")):
            problems.append(
                f"material.exe: PENELOPE database not found at {self.pdfiles}/ "
                "(material.exe must sit in pendbase/, next to pdfiles/)"
            )
        return problems

    def problems(self):
        """Human-readable list of everything that would stop a run."""
        return self.peneasy_problems() + self.material_problems()

    def require_peneasy(self):
        problems = self.peneasy_problems()
        if problems:
            raise ToolchainError(problems, ENV_PENEASY)
        return self.peneasy

    def require_material(self):
        problems = self.material_problems()
        if problems:
            raise ToolchainError(problems, ENV_MATERIAL)
        return self.material


class ToolchainError(RuntimeError):
    def __init__(self, problems, env_var):
        super().__init__(
            "\n  ".join(["PENELOPE toolchain problem:"] + list(problems))
            + f"\n  (pass the path explicitly or set ${env_var})"
        )
        self.problems = problems
