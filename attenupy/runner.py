"""Run penEasy over simulation directories.

A runner is any object with ``run(sim_dirs) -> RunSummary``.
:class:`LocalRunner` keeps a bounded pool of ``penEasy.x`` processes on the
current machine.

Coordination happens through files inside each simulation directory, so any
number of runners (e.g. several jobs on a shared filesystem) can work through the
same tree at once:

* ``.attenupy/claim`` -- created atomically before launching; whoever creates
  it owns the run. Claims left by dead processes on the same host are reclaimed
  automatically.
* ``.attenupy/result.json`` -- exit code and timings once the process ends.
* completion is judged by penEasy's own sign-off in ``penEasy.out``.

States (see :func:`sim_state`): ``pending``, ``running``, ``done``, ``failed``,
``unprepared``.
"""

import json
import logging
import os
import re
import signal
import socket
import subprocess
import time
from dataclasses import dataclass, field

from .results import is_finished

log = logging.getLogger(__name__)

STATE_DIR = ".attenupy"
PENDING, RUNNING, DONE, FAILED, UNPREPARED = "pending", "running", "done", "failed", "unprepared"


def _state_path(sim_dir, name):
    return os.path.join(sim_dir, STATE_DIR, name)


def _read_json(path):
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def _pid_alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _claim_is_stale(claim):
    return (claim is not None and claim.get("host") == socket.gethostname()
            and not _pid_alive(int(claim.get("pid", -1))))


# Input files named in penEasy.in: .spc/.geo lines and the material table.
_DECK_FILE_RE = re.compile(r"^\s+(?:\d+\s+)?(\S+\.(?:mat|geo|spc))\s", re.M)


def has_content(path):
    """True for a non-empty file. penEasy creates an *empty* file when an input
    such as a material is missing, so mere existence proves nothing."""
    try:
        return os.path.getsize(path) > 0
    except OSError:
        return False


def missing_inputs(sim_dir):
    """Files referenced by ``penEasy.in`` that are absent (or empty) in ``sim_dir``."""
    try:
        with open(os.path.join(sim_dir, "penEasy.in")) as f:
            deck = f.read()
    except OSError:
        return ["penEasy.in"]
    names = dict.fromkeys(_DECK_FILE_RE.findall(deck))
    return [n for n in names if not has_content(os.path.join(sim_dir, n))]


def sim_state(sim_dir):
    """Current state of one simulation directory."""
    if missing_inputs(sim_dir):
        return UNPREPARED
    if is_finished(sim_dir):
        return DONE
    claim_path = _state_path(sim_dir, "claim")
    if os.path.exists(claim_path):
        return PENDING if _claim_is_stale(_read_json(claim_path)) else RUNNING
    if os.path.exists(_state_path(sim_dir, "result.json")):
        return FAILED
    return PENDING


def reset(sim_dir):
    """Forget claims/results so the simulation is ``pending`` again."""
    for name in ("claim", "result.json"):
        try:
            os.remove(_state_path(sim_dir, name))
        except FileNotFoundError:
            pass


def _try_claim(sim_dir):
    os.makedirs(os.path.join(sim_dir, STATE_DIR), exist_ok=True)
    path = _state_path(sim_dir, "claim")
    if os.path.exists(path) and _claim_is_stale(_read_json(path)):
        os.remove(path)
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        return False
    with os.fdopen(fd, "w") as f:
        json.dump({"host": socket.gethostname(), "pid": os.getpid(),
                   "started": time.time()}, f)
    return True


def _raise_interrupt(signum, frame):
    raise KeyboardInterrupt


def _install_sigterm_handler():
    """Turn SIGTERM (e.g. a scheduler's time limit) into KeyboardInterrupt so running
    simulations are stopped and their claims released. Returns the previous
    handler, or None when not in the main thread."""
    try:
        previous = signal.getsignal(signal.SIGTERM)
        signal.signal(signal.SIGTERM, _raise_interrupt)
        return previous
    except ValueError:
        return None


@dataclass
class RunSummary:
    launched: int = 0
    succeeded: int = 0
    failed: list = field(default_factory=list)
    skipped: int = 0


def select(sim_dirs, *, retry_failed=False):
    """The directories a runner should attempt, in the given order."""
    wanted = {PENDING, FAILED} if retry_failed else {PENDING}
    return [d for d in sim_dirs if sim_state(d) in wanted]


class LocalRunner:
    """Run up to ``max_concurrent`` penEasy processes on this machine.

    ``executable`` is the penEasy.x to use; when ``None`` each directory's own
    ``penEasy.x`` is used instead (a folder that carries its own copy).
    """

    def __init__(self, executable=None, *, max_concurrent=None, limit=None,
                 retry_failed=False, poll_seconds=2.0):
        self.executable = executable
        self.max_concurrent = max_concurrent or os.cpu_count() or 1
        self.limit = limit
        self.retry_failed = retry_failed
        self.poll_seconds = poll_seconds

    def _exe_for(self, sim_dir):
        exe = self.executable or os.path.join(sim_dir, "penEasy.x")
        if not os.path.isfile(exe):
            raise FileNotFoundError(
                f"No penEasy.x for {sim_dir}: pass the executable (or set $ATTENUPY_PENEASY)")
        return exe

    def _launch(self, sim_dir):
        exe = self._exe_for(sim_dir)
        with open(os.path.join(sim_dir, "penEasy.in")) as fin, \
                open(os.path.join(sim_dir, "penEasy.out"), "w") as fout:
            return subprocess.Popen([exe], stdin=fin, stdout=fout,
                                    stderr=subprocess.STDOUT, cwd=sim_dir)

    @staticmethod
    def _finish(sim_dir, returncode, started, summary):
        with open(_state_path(sim_dir, "result.json"), "w") as f:
            json.dump({"returncode": returncode, "started": started,
                       "finished": time.time(), "host": socket.gethostname()}, f)
        os.remove(_state_path(sim_dir, "claim"))
        ok = returncode == 0 and is_finished(sim_dir)
        if ok:
            summary.succeeded += 1
            os.remove(_state_path(sim_dir, "result.json"))
        else:
            summary.failed.append(sim_dir)
        (log.info if ok else log.warning)("  %s (exit %s, %.0fs): %s", "done  " if ok else "FAILED",
                                          returncode, time.time() - started, sim_dir)

    def run(self, sim_dirs):
        todo = select(sim_dirs, retry_failed=self.retry_failed)
        summary = RunSummary()
        log.info("%d simulations to run, up to %d at a time", len(todo), self.max_concurrent)
        active = {}  # sim_dir -> (proc, started)

        def reap():
            for d, (proc, started) in list(active.items()):
                rc = proc.poll()
                if rc is not None:
                    del active[d]
                    self._finish(d, rc, started, summary)

        previous = _install_sigterm_handler()
        try:
            for sim_dir in todo:
                if self.limit is not None and summary.launched >= self.limit:
                    break
                while len(active) >= self.max_concurrent:
                    time.sleep(self.poll_seconds)
                    reap()
                if not _try_claim(sim_dir):
                    summary.skipped += 1  # another runner got there first
                    continue
                # The list was selected up front: another runner may have
                # finished (or failed) this one since. Re-check under the claim.
                result = _state_path(sim_dir, "result.json")
                if is_finished(sim_dir) or (os.path.exists(result) and not self.retry_failed):
                    os.remove(_state_path(sim_dir, "claim"))
                    summary.skipped += 1
                    continue
                if os.path.exists(result):
                    os.remove(result)  # retrying a failure
                try:
                    proc = self._launch(sim_dir)
                except Exception:
                    os.remove(_state_path(sim_dir, "claim"))
                    raise
                active[sim_dir] = (proc, time.time())
                summary.launched += 1
                log.info("  start (pid %d): %s", active[sim_dir][0].pid, sim_dir)
            while active:
                time.sleep(self.poll_seconds)
                reap()
        except KeyboardInterrupt:
            log.warning("Interrupted: stopping running simulations and releasing their claims")
            for d, (proc, _) in active.items():
                proc.terminate()
            for d, (proc, _) in active.items():
                proc.wait()
                reset(d)
            raise
        finally:
            if previous is not None:
                signal.signal(signal.SIGTERM, previous)
        log.info("Launched %d: %d done, %d failed, %d taken by other runners.",
                 summary.launched, summary.succeeded, len(summary.failed), summary.skipped)
        return summary
