"""Command line: ``attenupy <command> ...`` (or ``python -m attenupy``).

    attenupy doctor                                  check penEasy.x / material.exe
    attenupy prepare ROOT --design design.toml       create folders (+ materials)
    attenupy prepare ROOT --design d.toml --dry-run  just show what would be made
    attenupy materials ROOT                          build missing mat.mat / Aire.mat
    attenupy references ROOT                         add air runs for I0
    attenupy status ROOT
    attenupy run ROOT --max-concurrent 10
    attenupy collect ROOT -o results.csv
    attenupy reset ROOT --state failed

Executables come from ``--peneasy`` / ``--material`` or the environment
variables ``ATTENUPY_PENEASY`` / ``ATTENUPY_MATERIAL``.
"""

import argparse
import os
import sys
from collections import Counter

from . import log_to_console
from .design import load_design
from .naming import CompactNamer, DescriptiveNamer
from .results import DEFAULT_EXTRACTORS, SpecColumns, material_properties, EnergyRatio, path_column
from .runner import FAILED, RUNNING, LocalRunner, reset, sim_state
from .study import Study
from .toolchain import Toolchain

_NAMERS = {"compact": CompactNamer, "descriptive": DescriptiveNamer}


def _study(args):
    return Study(args.root, toolchain=Toolchain.from_env(args.peneasy, args.material),
                 namer=_NAMERS[getattr(args, "namer", "compact")]())


def _cmd_doctor(args):
    tc = Toolchain.from_env(args.peneasy, args.material)
    print(f"penEasy.x    : {tc.peneasy or '(not set)'}")
    print(f"material.exe : {tc.material or '(not set)'}")
    problems = tc.problems()
    for p in problems:
        print(f"  ✗ {p}")
    if not problems:
        print("  ✓ both executables found and runnable on this machine")
    return 1 if problems else 0


def _cmd_prepare(args):
    study = _study(args)
    specs = list(load_design(args.design, constituents=study.constituents, beams=study.beams))
    plan = study.plan(specs)
    by_beam = Counter(study.beams[s.beam].folder for s, _ in plan)
    by_size = Counter(s.mixture.n_constituents for s, _ in plan if s.mixture)
    print(f"{len(plan)} simulations: " + ", ".join(f"{k}: {v}" for k, v in sorted(by_beam.items())))
    print("constituents per mixture: " + ", ".join(f"{k}: {v}" for k, v in sorted(by_size.items())))
    if args.dry_run:
        for _, path in plan[:10]:
            print("  " + os.path.relpath(path, study.root))
        if len(plan) > 10:
            print(f"  ... ({len(plan) - 10} more)")
        return 0
    study.prepare(specs, build_materials=not args.no_materials, workers=args.workers)
    return 0


def _cmd_materials(args):
    _study(args).build_materials(overwrite=args.overwrite, workers=args.workers)
    return 0


def _cmd_references(args):
    _study(args).prepare_references(args.beam or None)
    return 0


def _cmd_status(args):
    study = _study(args)
    table = Counter()
    for sim in study.simulations(include_references=True):
        table[(os.path.relpath(os.path.dirname(sim.path), study.root), sim_state(sim.path))] += 1
    states = ["unprepared", "pending", "running", "done", "failed"]
    groups = sorted({g for g, _ in table})
    print(f"{'folder':<12}" + "".join(f"{s:>11}" for s in states))
    for g in groups:
        print(f"{g:<12}" + "".join(f"{table[(g, s)]:>11}" for s in states))
    if any(s == "unprepared" for _, s in table):
        print("unprepared = inputs missing (e.g. materials): run `attenupy materials ROOT`")
    return 0


def _cmd_run(args):
    study = _study(args)
    runner = LocalRunner(study.toolchain.peneasy, max_concurrent=args.max_concurrent,
                         limit=args.limit, retry_failed=args.retry_failed)
    summary = study.run(runner, beams=args.beam)
    return 1 if summary.failed else 0


def _cmd_collect(args):
    extractors = ([SpecColumns(), material_properties, EnergyRatio()] if args.tidy
                  else list(DEFAULT_EXTRACTORS))
    if args.with_path:
        extractors.append(path_column)
    _study(args).collect(extractors, output=args.output,
                         only_finished=not args.include_unfinished)
    return 0


def _cmd_reset(args):
    n = 0
    for sim in _study(args).simulations(include_references=True):
        if sim_state(sim.path) in args.state:
            reset(sim.path)
            n += 1
    print(f"Reset {n} simulation(s) to pending.")
    return 0


def build_parser():
    p = argparse.ArgumentParser(prog="attenupy", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--peneasy", help="path to penEasy.x (or $ATTENUPY_PENEASY)")
    p.add_argument("--material", help="path to material.exe (or $ATTENUPY_MATERIAL)")
    sub = p.add_subparsers(dest="command", required=True)

    def cmd(name, func, help_, root=True):
        sp = sub.add_parser(name, help=help_)
        if root:
            sp.add_argument("root", help="study directory")
        sp.set_defaults(func=func)
        return sp

    cmd("doctor", _cmd_doctor, "check the PENELOPE executables", root=False)

    sp = cmd("prepare", _cmd_prepare, "create simulation folders from a design file")
    sp.add_argument("--design", required=True, help="design file (.json / .toml)")
    sp.add_argument("--namer", choices=sorted(_NAMERS), default="compact")
    sp.add_argument("--no-materials", action="store_true",
                    help="write inputs only; run `attenupy materials` later")
    sp.add_argument("--dry-run", action="store_true")
    sp.add_argument("--workers", type=int, default=None, help="parallel material.exe runs")

    sp = cmd("materials", _cmd_materials, "build missing material files")
    sp.add_argument("--overwrite", action="store_true")
    sp.add_argument("--workers", type=int, default=None)

    sp = cmd("references", _cmd_references, "add air reference runs (I0)")
    sp.add_argument("--beam", action="append", help="beam key/alias (default: all used)")

    cmd("status", _cmd_status, "count pending/running/done/failed")

    sp = cmd("run", _cmd_run, "run pending simulations")
    sp.add_argument("--max-concurrent", type=int, default=None)
    sp.add_argument("--limit", type=int, default=None, help="launch at most N")
    sp.add_argument("--beam", action="append", help="only these beams")
    sp.add_argument("--retry-failed", action="store_true")

    sp = cmd("collect", _cmd_collect, "harvest results into a CSV")
    sp.add_argument("-o", "--output", default="results.csv")
    sp.add_argument("--tidy", action="store_true",
                    help="one w_<constituent> column each instead of molecule_i/ratio_i")
    sp.add_argument("--with-path", action="store_true")
    sp.add_argument("--include-unfinished", action="store_true")

    sp = cmd("reset", _cmd_reset, "mark failed/stuck simulations pending again")
    sp.add_argument("--state", action="append", choices=[FAILED, RUNNING], required=True)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    log_to_console("INFO")
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()
