"""
Process-local tracking of the Slurm jobs and local driver processes that one
`run_brb` invocation has started, so a fatal crash can cancel exactly its own
work. Deliberately in-memory and per-process: `run_brb -s illumina` and
`run_brb -s aviti` can run side by side on one host, and neither may cancel
the other's jobs.
"""

import json
import os
import re
import threading
import time
from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path

from BRB.logger import log


class GroupHandle:
    """One library-group's in-flight state: its Slurm job IDs and drivers."""

    def __init__(self, outputDir):
        self.outputDir = str(outputDir)
        self._lock = threading.Lock()
        self._job_ids = []
        self._processes = []

    def add_job_ids(self, ids):
        """Record newly-seen Slurm job IDs. Returns only the ones that were new."""
        with self._lock:
            new = [i for i in ids if i not in self._job_ids]
            self._job_ids.extend(new)
        if new:
            log.debug(f"Tracking Slurm job(s) {','.join(new)} for {self.outputDir}")
        return new

    def add_process(self, proc):
        with self._lock:
            self._processes.append(proc)

    def remove_process(self, proc):
        with self._lock:
            if proc in self._processes:
                self._processes.remove(proc)

    def snapshot(self):
        """(job_ids, processes) copies, safe to iterate outside the lock."""
        with self._lock:
            return list(self._job_ids), list(self._processes)


class JobRegistry:
    """
    All groups currently in flight for one flowcell, in one `run_brb` process.
    Construct it in `runFlowcell` and pass it down; never a module-level global.
    """

    def __init__(self):
        self._lock = threading.Lock()
        self._groups = {}
        self._aborted = False

    @property
    def aborted(self):
        return self._aborted

    def abort(self):
        """Latch 'this flowcell is being torn down'; workers must stop dispatching."""
        self._aborted = True

    def register_group(self, outputDir):
        key = str(outputDir)
        with self._lock:
            handle = self._groups.get(key)
            if handle is None:
                handle = GroupHandle(key)
                self._groups[key] = handle
            return handle

    def unregister_group(self, outputDir):
        with self._lock:
            return self._groups.pop(str(outputDir), None)

    def active_groups(self):
        with self._lock:
            return list(self._groups.values())


_currentHandle = ContextVar("BRB_group_handle", default=None)


def currentHandle():
    """The GroupHandle bound to the calling thread, or None."""
    return _currentHandle.get()


@contextmanager
def bindGroup(handle):
    """
    Bind `handle` for the duration of the block, so the pipeline functions'
    `runManagedSubprocess` calls can find it without threading a parameter
    through ten unrelated signatures. Scoped to the calling thread, so it is
    not shared mutable state between pool workers.
    """
    token = _currentHandle.set(handle)
    try:
        yield handle
    finally:
        _currentHandle.reset(token)


MARKER_NAME = "running.pid"

MARKER_ABSENT = "absent"
MARKER_LIVE = "live"
MARKER_ABANDONED = "abandoned"
MARKER_CORRUPT = "corrupt"


def markerPath(outputDir):
    return Path(outputDir) / MARKER_NAME


def _atomicWrite(outputDir, payload):
    """
    Write the marker via a temp file + os.replace, so a crash mid-write can
    never leave a half-written marker behind (see markerState's CORRUPT case).
    """
    target = markerPath(outputDir)
    tmp = target.parent / (MARKER_NAME + ".tmp")
    tmp.write_text(json.dumps(payload))
    os.replace(tmp, target)
    return payload


def writeMarker(outputDir, pid=None, job_ids=None, cancelled=False):
    """Claim ownership of `outputDir` for this process."""
    return _atomicWrite(
        outputDir,
        {
            "pid": os.getpid() if pid is None else pid,
            "started": time.time(),
            "job_ids": list(job_ids or []),
            "cancelled": bool(cancelled),
        },
    )


def readMarker(outputDir):
    """The parsed marker dict, or None if it is absent or not a JSON object."""
    try:
        data = json.loads(markerPath(outputDir).read_text())
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def updateMarkerJobIds(outputDir, job_ids):
    """Rewrite the marker's job_ids, keeping pid/started/cancelled as they were."""
    existing = readMarker(outputDir) or {}
    return _atomicWrite(
        outputDir,
        {
            "pid": existing.get("pid", os.getpid()),
            "started": existing.get("started", time.time()),
            "job_ids": list(job_ids),
            "cancelled": bool(existing.get("cancelled", False)),
        },
    )


def markMarkerCancelled(outputDir):
    """
    Flag the marker as cancelled *before* the scancel/SIGTERM goes out, so a
    crash handler that dies partway still leaves a readable marker that says
    what happened.
    """
    existing = readMarker(outputDir)
    if existing is None:
        return None
    existing["cancelled"] = True
    return _atomicWrite(outputDir, existing)


def clearMarker(outputDir):
    markerPath(outputDir).unlink(missing_ok=True)


def _pidAlive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        # Someone else's process: alive, just not ours to signal.
        return True
    except OSError:
        return True
    return True


def markerState(outputDir):
    """
    Classify `outputDir`'s ownership marker as (state, marker).

    ABSENT     -> nobody owns it, dispatch normally.
    LIVE       -> another live process owns it, skip this group this pass.
    ABANDONED  -> the owner is gone, remove the marker and dispatch.
    CORRUPT    -> unreadable; it cannot prove the owner is dead, so skip the
                  group and let an operator clear it by hand.
    """
    path = markerPath(outputDir)
    if not path.exists():
        return (MARKER_ABSENT, None)
    marker = readMarker(outputDir)
    if marker is None or not isinstance(marker.get("pid"), int):
        return (MARKER_CORRUPT, None)
    if _pidAlive(marker["pid"]):
        return (MARKER_LIVE, marker)
    return (MARKER_ABANDONED, marker)


# snakemake's cluster executors log one of:
#   Submitted job 12 with external jobid 'Submitted batch job 1234567'.
#   Submitted job 12 with external jobid '1234567'.        (sbatch --parsable)
# The quoted payload is whatever the submit command printed, so take the
# first run of digits in it as the Slurm job ID.
_EXTERNAL_JOBID_RE = re.compile(r"external jobid[:\s]*'([^']*)'")
_FIRST_INT_RE = re.compile(r"\d+")


def parseJobIds(line):
    """Slurm job IDs announced on one line of snakemake driver output."""
    found = []
    for payload in _EXTERNAL_JOBID_RE.findall(line):
        match = _FIRST_INT_RE.search(payload)
        if match:
            found.append(match.group(0))
    return found
