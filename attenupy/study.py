"""A :class:`Study` ties the pieces together for one directory tree.

    study = Study("runs/pe_composites", toolchain=Toolchain.from_env())
    study.prepare(RandomSampler(...))     # folders, spec.json, penEasy.in, .geo, .spc, materials
    study.prepare_references()            # air runs to measure I0 (only for new beams)
    study.run(max_concurrent=10)          # penEasy.x, resumable, multi-job safe
    rows = study.collect(output="results.csv")

Layout: ``<root>/<beam.folder>/<namer.name(spec)>/``, plus ``<root>/.attenupy/``
for the material cache and provenance. Every component is swappable through
the constructor: ``geometry``, ``deck``, ``namer``, ``materials`` (builder),
``constituents`` and ``beams`` registries.
"""

import hashlib
import json
import logging
import os
import shutil
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

from . import __version__
from .beams import DEFAULT_BEAMS
from .constituents import DEFAULT_CONSTITUENTS
from .geometry import SlabGeometry
from .materials import PenelopeMaterials
from .naming import CompactNamer
from .results import (DEFAULT_EXTRACTORS, CollectContext, Simulation, collect,
                      find_simulations, is_finished, load_spec,
                      read_energy_deposition, write_csv, DETECTOR_MATERIAL,
                      ENERGY_TALLY)
from .runner import LocalRunner, has_content, sim_state
from .spec import SimulationSpec
from .template import PeneasyDeck
from .toolchain import Toolchain

log = logging.getLogger(__name__)

REFERENCE_DIRNAME = "_reference"
SLAB_MAT, AIR_MAT = "mat.mat", "Aire.mat"


class Study:
    def __init__(self, root, *, toolchain=None, geometry=None, deck=None, namer=None,
                 materials=None, constituents=None, beams=None):
        self.root = os.path.abspath(os.path.expanduser(root))
        self.toolchain = toolchain or Toolchain.from_env()
        self.geometry = geometry or SlabGeometry()
        self.deck = deck or PeneasyDeck()
        self.namer = namer or CompactNamer()
        self.constituents = constituents or DEFAULT_CONSTITUENTS
        self.beams = beams or DEFAULT_BEAMS
        self._materials = materials

    # ------------------------------------------------------------------
    # Components
    # ------------------------------------------------------------------
    @property
    def materials(self):
        """The material builder (default: material.exe from the toolchain)."""
        if self._materials is None:
            self._materials = PenelopeMaterials(self.toolchain.require_material())
        return self._materials

    @property
    def state_dir(self):
        return os.path.join(self.root, ".attenupy")

    def path_for(self, spec):
        beam = self.beams[spec.beam]
        name = REFERENCE_DIRNAME if spec.is_reference else self.namer.name(spec)
        return os.path.join(self.root, beam.folder, name)

    # ------------------------------------------------------------------
    # Planning: validate everything before writing anything
    # ------------------------------------------------------------------
    def _validate(self, spec):
        beam = self.beams[spec.beam]
        spec = spec.replace(beam=beam.key)
        if hasattr(beam.spectrum, "validate"):
            beam.spectrum.validate()
        if spec.mixture is not None:
            for name, _ in spec.mixture.components():
                self.constituents[name]
        if hasattr(self.geometry, "validate"):
            self.geometry.validate(spec)
        for fname in (beam.spectrum_filename, self.geometry.filename(spec)):
            if len(fname) > 80:
                raise ValueError(f"penEasy limits file names to 80 chars: {fname!r}")
        return spec

    def plan(self, specs):
        """``[(spec, path), ...]`` for ``specs``, validated and de-duplicated.

        Identical specs are kept once. Different specs that the namer maps to
        the same folder get ``__2``, ``__3``... suffixes, as do specs whose
        folder already exists on disk holding a different simulation.
        """
        angle = getattr(self.deck, "beam_half_angle_deg", None)
        if angle is not None and hasattr(self.geometry, "check_beam"):
            self.geometry.check_beam(float(angle))
        planned, problems, taken, duplicates = [], [], {}, 0
        for i, spec in enumerate(specs):
            try:
                spec = self._validate(spec)
                base = self.path_for(spec)
            except (KeyError, ValueError, FileNotFoundError) as e:
                problems.append(f"spec #{i} ({spec.beam}, {spec.thickness_mm} mm): "
                                f"{e.args[0] if e.args else e}")
                continue
            path, n = base, 1
            while True:
                existing = taken.get(path)
                if existing is None and os.path.isdir(path):
                    try:
                        existing = load_spec(path).replace(beam=spec.beam)
                    except (ValueError, KeyError, OSError):
                        existing = "unknown"
                if existing is None or existing == spec:
                    break
                n += 1
                path = f"{base}__{n}"
            if taken.get(path) == spec:
                duplicates += 1
                continue
            taken[path] = spec
            planned.append((spec, path))
        if problems:
            more = f"\n  ... and {len(problems) - 10} more" if len(problems) > 10 else ""
            raise ValueError("Invalid design:\n  " + "\n  ".join(problems[:10]) + more)
        if duplicates:
            log.info("Dropped %d duplicate spec(s).", duplicates)
        return planned

    # ------------------------------------------------------------------
    # Preparing directories
    # ------------------------------------------------------------------
    def prepare(self, specs, *, build_materials=True, overwrite=False,
                copy_executable=False, workers=None):
        """Create simulation directories for ``specs``. Returns their paths.

        With ``build_materials=False`` only the text inputs are written (e.g.
        on a laptop); run :meth:`build_materials` later where material.exe is.
        """
        plan = self.plan(specs)
        if build_materials:
            self.materials  # fail now, not after writing half the tree
        if copy_executable:
            self.toolchain.require_peneasy()
        created = 0
        for spec, path in plan:
            if overwrite or not os.path.isfile(os.path.join(path, "penEasy.in")):
                self._write_inputs(spec, path, copy_executable)
                created += 1
        self._write_provenance()
        log.info("Prepared %d simulation folder(s) under %s (%d already existed).",
                 created, self.root, len(plan) - created)
        if build_materials:
            self._build_materials(plan, overwrite=overwrite, workers=workers)
        return [path for _, path in plan]

    def _write_inputs(self, spec, path, copy_executable=False):
        os.makedirs(path, exist_ok=True)
        beam = self.beams[spec.beam]
        spec.save(path)
        beam.spectrum.write(os.path.join(path, beam.spectrum_filename))
        geo = self.geometry.filename(spec)
        self.geometry.write(spec, os.path.join(path, geo))
        self.deck.write(path, beam.spectrum_filename, geo, seed=spec.seed,
                        slab_material=AIR_MAT if spec.is_reference else SLAB_MAT,
                        beam=beam.key, kvp=beam.kvp, thickness_mm=spec.thickness_mm)
        if copy_executable:
            dst = os.path.join(path, "penEasy.x")
            shutil.copyfile(self.toolchain.peneasy, dst)
            os.chmod(dst, 0o755)

    def _write_provenance(self):
        os.makedirs(self.state_dir, exist_ok=True)
        deck = {k: v for k, v in vars(self.deck).items() if k != "template"}
        deck["template_sha1"] = hashlib.sha1(self.deck.template.encode()).hexdigest()
        info = {
            "attenupy_version": __version__,
            "geometry": repr(self.geometry),
            "deck": deck,
            "namer": type(self.namer).__name__,
            "toolchain": {"peneasy": self.toolchain.peneasy,
                          "material": self.toolchain.material},
        }
        with open(os.path.join(self.state_dir, "study.json"), "w") as f:
            json.dump(info, f, indent=2)

    # ------------------------------------------------------------------
    # Materials
    # ------------------------------------------------------------------
    def _material_key(self, mixture):
        comps = {n: [self.constituents[n].formula, self.constituents[n].density]
                 for n, _ in mixture.components()}
        blob = json.dumps({"mixture": mixture.to_dict(), "constituents": comps}, sort_keys=True)
        return hashlib.sha1(blob.encode()).hexdigest()[:16]

    def _cached(self, name, build):
        """Path of a cached material, building it first if needed."""
        cache = os.path.join(self.state_dir, "materials")
        os.makedirs(cache, exist_ok=True)
        path = os.path.join(cache, name)
        if not os.path.isfile(path):
            tmp = f"{path}.tmp{os.getpid()}"
            build(tmp)
            os.replace(tmp, path)
        return path

    def _build_materials(self, plan, *, overwrite=False, workers=None):
        builder = self.materials
        todo = [(spec, path) for spec, path in plan
                if overwrite or not has_content(os.path.join(path, AIR_MAT))
                or (not spec.is_reference and not has_content(os.path.join(path, SLAB_MAT)))]
        if not todo:
            return
        air = self._cached(AIR_MAT, builder.build_air)
        unique = {self._material_key(s.mixture): s.mixture for s, _ in todo if not s.is_reference}
        missing = sum(not os.path.isfile(os.path.join(self.state_dir, "materials", f"{k}.mat"))
                      for k in unique)
        log.info("Materials for %d folder(s): %d unique, %d to build with material.exe...",
                 len(todo), len(unique), missing)

        def build_one(item):
            key, mixture = item
            return key, self._cached(
                f"{key}.mat", lambda dst: builder.build(mixture, dst, self.constituents))

        with ThreadPoolExecutor(max_workers=workers or os.cpu_count() or 1) as pool:
            built = dict(pool.map(build_one, unique.items()))
        for spec, path in todo:
            shutil.copyfile(air, os.path.join(path, AIR_MAT))
            if not spec.is_reference:
                shutil.copyfile(built[self._material_key(spec.mixture)],
                                os.path.join(path, SLAB_MAT))
        log.info("Materials ready.")

    def build_materials(self, *, overwrite=False, workers=None):
        """Build missing ``mat.mat`` / ``Aire.mat`` for every folder in the study.

        Materials come from the cache in ``.attenupy/materials`` (keyed by
        composition and constituent data); ``overwrite=True`` re-copies them.
        Delete that directory to force material.exe to run again.
        """
        plan = [(s.spec.replace(beam=s.beam.key if s.beam else s.spec.beam), s.path)
                for s in self.simulations(include_references=True)]
        self._build_materials(plan, overwrite=overwrite, workers=workers)

    # ------------------------------------------------------------------
    # References (I0)
    # ------------------------------------------------------------------
    def prepare_references(self, beams=None, *, thickness_mm=0.30, build_materials=True):
        """Create one air run per beam (default: every beam used in the study).

        Its detector energy is the I0 that :meth:`collect` divides by; it takes
        precedence over a beam's configured ``reference_energy``.
        """
        if beams is None:
            beams = sorted({s.spec.beam for s in self.simulations()})
        specs = [SimulationSpec.reference(self.beams[b].key, thickness_mm) for b in beams]
        return self.prepare(specs, build_materials=build_materials)

    def measured_references(self):
        """``{beam_key: I0}`` from finished reference runs."""
        out = {}
        for sim in self.simulations(include_references=True):
            if sim.spec.is_reference and is_finished(sim.path):
                tally = read_energy_deposition(sim.file(ENERGY_TALLY))
                if tally and DETECTOR_MATERIAL in tally:
                    out[sim.beam.key if sim.beam else sim.spec.beam] = tally[DETECTOR_MATERIAL][0]
        return out

    # ------------------------------------------------------------------
    # Inspecting, running, collecting
    # ------------------------------------------------------------------
    def simulations(self, include_references=False):
        if not os.path.isdir(self.root):
            return []
        return find_simulations(self.root, self.beams, include_references)

    def status(self):
        """``Counter`` of simulation states (pending/running/done/failed)."""
        return Counter(sim_state(s.path) for s in self.simulations(include_references=True))

    def run(self, runner=None, *, beams=None, **runner_kwargs):
        """Run every pending simulation (reference runs first).

        ``runner`` is any object with ``run(dirs)``; by default a
        :class:`LocalRunner` built from ``runner_kwargs``. ``beams`` restricts
        the run to some beams (keys or aliases).
        """
        if runner is None:
            runner = LocalRunner(self.toolchain.peneasy, **runner_kwargs)
        sims = self.simulations(include_references=True)
        if beams:
            wanted = {self.beams[b].key for b in beams}
            sims = [s for s in sims if s.beam and s.beam.key in wanted]
        refs = [s.path for s in sims if s.spec.is_reference]
        return runner.run(refs + [s.path for s in sims if not s.spec.is_reference])

    def collect(self, extractors=None, output=None, *, only_finished=True):
        """Rows for finished simulations; also written to ``output`` if given."""
        ctx = CollectContext(self.beams, self.measured_references(), self.constituents)
        rows = collect(self.simulations(), extractors or DEFAULT_EXTRACTORS, ctx,
                       only_finished=only_finished)
        if output:
            write_csv(rows, output)
            log.info("Wrote %s", output)
        return rows


__all__ = ["Study", "Simulation"]
