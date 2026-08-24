"""
Process-local tracking of the Slurm jobs and local driver processes that one
`run_brb` invocation has started, so a fatal crash can cancel exactly its own
work. Deliberately in-memory and per-process: `run_brb -s illumina` and
`run_brb -s aviti` can run side by side on one host, and neither may cancel
the other's jobs.
"""

import threading
from contextlib import contextmanager
from contextvars import ContextVar

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
