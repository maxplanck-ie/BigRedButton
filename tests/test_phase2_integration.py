"""
Phase 2 integration: a crash mid-flowcell must cancel every in-flight group's
Slurm jobs, and the next pass must redo exactly the groups that were not
finished.
"""

import configparser
import os
import threading
import time

import pytest

from BRB import PushButton, jobtrack


@pytest.fixture
def flowcell(tmp_path, monkeypatch):
    """
    Three groups in two projects, each with its own outputDir. Pipeline
    functions are stubbed: they submit a fake Slurm job, honour analysis.done
    the way the real ones do, and record every dispatch.
    """
    config = configparser.ConfigParser()
    config["Paths"] = {
        "baseData": str(tmp_path / "base"),
        "groupData": str(tmp_path / "data"),
        "seqFacDir": str(tmp_path / "seqfac"),
        "bioinfoCoreDir": str(tmp_path / "bioinfocore"),
    }
    config["Options"] = {
        "runID": "20260824_AV999999_1234567_lanes_1_2",
        "sequencerType": "Aviti",
        "snakemakeWorkflowBaseDir": "/snakepipes",
    }
    config["external"] = {"LibraryTypes": "", "LibraryProtocols": ""}

    dirs = {}
    for name in ("alpha", "beta", "gamma"):
        d = tmp_path / "out" / name
        d.mkdir(parents=True)
        dirs[name] = d

    items = [
        PushButton.WorkItem(
            project="1_Doe_Smith",
            group="smith",
            pipeline="RNA",
            organism=("mouse", "GRCm38", "/yaml/GRCm38.yaml"),
            libraryType="stranded mRNA-Seq",
            tuples=[["18L001", "alpha", "TruSeq", False]],
        ),
        PushButton.WorkItem(
            project="1_Doe_Smith",
            group="smith",
            pipeline="RNA",
            organism=("human", "GRCh38", "/yaml/GRCh38.yaml"),
            libraryType="ChIP-Seq",
            tuples=[["18L002", "beta", "TruSeq", False]],
        ),
        PushButton.WorkItem(
            project="2_Roe_Jones",
            group="jones",
            pipeline="RNA",
            organism=("mouse", "GRCm38", "/yaml/GRCm38.yaml"),
            libraryType="stranded mRNA-Seq",
            tuples=[["18L003", "gamma", "TruSeq", False]],
        ),
    ]
    byName = {item.tuples[0][1]: dirs[item.tuples[0][1]] for item in items}

    monkeypatch.setattr(
        PushButton,
        "createPath",
        lambda config, group, project, org_label, libraryType, tuples: str(
            byName[tuples[0][1]]
        ),
    )
    monkeypatch.setattr(os.path, "exists", _realExistsExceptProjectDirs)
    monkeypatch.setattr(
        PushButton,
        "GetResults",
        lambda config, project, libraries: (
            [i for i in items if i.project == project],
            [],
        ),
    )
    monkeypatch.setattr(
        PushButton.BRB.ET,
        "phoneHome",
        lambda config, outputDir, pipeline, tuples, org, project, libType: [
            project,
            org,
            libType,
            pipeline,
            "OK",
            "updated",
        ],
    )
    parkour = {"1_Doe_Smith": {}, "2_Roe_Jones": {}}
    return config, parkour, dirs, items


def _realExistsExceptProjectDirs(path):
    # runFlowcell's "lives on another lane" guard: every project is present.
    if "/Project_" in str(path):
        return True
    return os.path.lexists(path)


def _stubPipeline(
    dirs, dispatched, behaviour, monkeypatch, registry=None, started=None
):
    """
    behaviour: {sampleName: "ok" | "crash" | "wait-for-abort"}.
    Emulates the real pipelines' analysis.done short-circuit and registers a
    fake Slurm job id on the bound handle.

    "wait-for-abort" exists to make the crash-cancellation tests
    deterministic: with a bounded thread pool, a group queued right after
    the crashing one can otherwise start racing against runFlowcell's own
    crash-handling (which fires from a different thread, see
    `cancelAllGroups`/`executor.shutdown(cancel_futures=True)`), so whether
    it finishes before or after teardown reaches it is timing-dependent. A
    "wait-for-abort" group registers its (fake) Slurm job id immediately --
    same as "ok" -- then blocks (bounded, polling) until `registry.aborted`
    is set, then raises instead of touching analysis.done. Registering
    first, synchronously, before any waiting happens, is what makes
    `cancelAllGroups` reliably find it (and its job id) still in
    `registry.active_groups()` and actually cancel it, standing in for a
    real Slurm driver that submitted its job and then got scancelled mid-run
    rather than one that raced to completion first. It must not write to the
    marker again *after* the wait: that would race `cancelAllGroups`, which
    clears the marker from the main thread as soon as it cancels this same
    group, and whichever write lands last would win.

    `started`: optional {sampleName: threading.Event}. If given, the matching
    event is set the instant that group's stub is entered -- before it does
    anything else, including the analysis.done short-circuit. This lets a
    test tell "the pool never got around to starting this group" (future
    genuinely cancelled by `cancel_futures`) apart from "it started and is
    still running/blocked", which `registry.active_groups()` cannot: a group
    submitted right after a crashing one can sit in the pool's queue for an
    unbounded (if usually tiny) amount of wall time before its worker thread
    executes a single line of `runOneGroup`, so polling any state the group
    itself would produce is inherently racy until it's known to have started.
    """

    def fake(config, group, project, organism, libraryType, tuples):
        name = tuples[0][1]
        if started is not None and name in started:
            started[name].set()
        outputDir = str(dirs[name])
        if os.path.lexists(os.path.join(outputDir, "analysis.done")):
            return outputDir, 0, False
        what = behaviour.get(name, "ok")
        dispatched.append(name)
        handle = jobtrack.currentHandle()
        if handle is not None:
            handle.add_job_ids([f"90{len(dispatched)}"])
            jobtrack.updateMarkerJobIds(outputDir, handle.snapshot()[0])
        if what == "crash":
            raise RuntimeError(f"driver for {name} died")
        if what == "wait-for-abort" and registry is not None:
            deadline = time.monotonic() + 5.0
            while not registry.aborted:
                if time.monotonic() > deadline:
                    raise AssertionError(
                        f"{name}: timed out waiting for registry.aborted"
                    )
                time.sleep(0.005)
            raise RuntimeError(f"driver for {name} killed by cancellation")
        PushButton.touchDone(outputDir)
        return outputDir, 0, False

    monkeypatch.setattr(PushButton, "RNA", fake)


def _waitFor(predicate, timeout=5.0, what="condition"):
    """Poll `predicate` until it's true, or raise after `timeout` seconds."""
    deadline = time.monotonic() + timeout
    while not predicate():
        if time.monotonic() > deadline:
            raise AssertionError(f"timed out waiting for {what}")
        time.sleep(0.005)


def _settleCrashedGroup(
    dirs, name, startedEvent, notStartedGrace=0.3, stabilize=0.1, timeout=5.0
):
    """
    After a crashing `runFlowcell` call, wait for a "wait-for-abort" group's
    fate to be fully settled on disk before inspecting it.

    Two outcomes are both correct, and which one happens is a genuine race
    between this group's worker thread and `executor.shutdown(wait=False,
    cancel_futures=True)` in runFlowcell (see `_stubPipeline`'s docstring):

    - The pool never started it (`cancel_futures` won the race): its future
      was cancelled before `runOneGroup` ever ran, so no marker was ever
      written for it. We wait `notStartedGrace` seconds for `startedEvent`
      as insurance against a slow-scheduled worker thread firing *after* we
      already concluded "never started" -- then confirm it's still unset.
    - It started (the worker thread won the race): both `cancelAllGroups`
      (main thread) and this group's own `finally: clearMarker(...)`
      (worker thread) clear the same marker, and either can run first --
      whichever hasn't fired yet when we first observe MARKER_ABSENT can
      still fire moments later (itself a no-op, but a caller that has
      already written something new to that path in the meantime would
      have it silently deleted from under it). So we don't trust the first
      MARKER_ABSENT reading: we wait for it to *stay* absent for a
      `stabilize` window before calling it settled.
    """
    if not startedEvent.wait(timeout=notStartedGrace):
        assert not startedEvent.is_set(), (
            f"{name}: worker started only after the not-started grace period"
        )
        return
    deadline = time.monotonic() + timeout
    while True:
        _waitFor(
            lambda: jobtrack.markerState(dirs[name])[0] == jobtrack.MARKER_ABSENT,
            timeout=max(deadline - time.monotonic(), 0.001),
            what=f"{name} marker to clear",
        )
        time.sleep(stabilize)
        if jobtrack.markerState(dirs[name])[0] == jobtrack.MARKER_ABSENT:
            return
        if time.monotonic() > deadline:
            raise AssertionError(f"{name}: marker never stayed settled")


class TestCrashAndRestart:
    def test_crash_cancels_all_groups_and_leaves_no_markers(
        self, flowcell, monkeypatch
    ):
        config, parkour, dirs, _items = flowcell
        dispatched = []
        scancelled = []
        monkeypatch.setattr(
            jobtrack.subprocess,
            "run",
            lambda argv, **kw: scancelled.extend(argv[1:]),
        )
        # Track exactly which job id(s) each group's handle picks up, keyed
        # by outputDir, so we can assert precisely on beta's own job id below
        # instead of merely "something got scancelled".
        jobIdsByGroup = {}
        origAddJobIds = jobtrack.GroupHandle.add_job_ids

        def trackingAddJobIds(self, ids):
            new = origAddJobIds(self, ids)
            if new:
                jobIdsByGroup.setdefault(self.outputDir, []).extend(new)
            return new

        monkeypatch.setattr(jobtrack.GroupHandle, "add_job_ids", trackingAddJobIds)
        monkeypatch.setattr(PushButton, "POOL_SIZE", 3)
        registry = jobtrack.JobRegistry()
        started = {"alpha": threading.Event(), "gamma": threading.Event()}
        _stubPipeline(
            dirs,
            dispatched,
            {"beta": "crash", "alpha": "wait-for-abort", "gamma": "wait-for-abort"},
            monkeypatch,
            registry=registry,
            started=started,
        )

        with pytest.raises(PushButton.GroupDispatchError) as excinfo:
            PushButton.runFlowcell(config, parkour, registry)
        assert "driver for beta died" in str(excinfo.value.__cause__)

        _settleCrashedGroup(dirs, "alpha", started["alpha"])
        _settleCrashedGroup(dirs, "gamma", started["gamma"])
        _waitFor(lambda: registry.active_groups() == [], what="registry to drain")
        assert registry.active_groups() == []
        for d in dirs.values():
            assert jobtrack.markerState(d)[0] == jobtrack.MARKER_ABSENT
        # beta is the group that actually crashed: after the Finding 3 fix it
        # stays registered through the crash, so cancelAllGroups must scancel
        # its own tracked job id specifically, not just some sibling's.
        betaJobIds = jobIdsByGroup[str(dirs["beta"])]
        assert betaJobIds
        for jobId in betaJobIds:
            assert jobId in scancelled

    def test_restart_redoes_only_unfinished_groups(self, flowcell, monkeypatch):
        config, parkour, dirs, _items = flowcell
        monkeypatch.setattr(jobtrack.subprocess, "run", lambda *a, **k: None)
        monkeypatch.setattr(PushButton, "POOL_SIZE", 1)

        firstPass = []
        registry1 = jobtrack.JobRegistry()
        gammaStarted = threading.Event()
        _stubPipeline(
            dirs,
            firstPass,
            {"beta": "crash", "gamma": "wait-for-abort"},
            monkeypatch,
            registry=registry1,
            started={"gamma": gammaStarted},
        )
        with pytest.raises(RuntimeError):
            PushButton.runFlowcell(config, parkour, registry1)
        # gamma is the "wait-for-abort" group -- settle its fate (never
        # started vs. started-then-cancelled) before treating the crash's
        # on-disk aftermath as final.
        _settleCrashedGroup(dirs, "gamma", gammaStarted)

        finished = {
            name
            for name, d in dirs.items()
            if os.path.lexists(os.path.join(str(d), "analysis.done"))
        }
        assert "beta" not in finished

        secondPass = []
        _stubPipeline(dirs, secondPass, {}, monkeypatch)
        msg = PushButton.runFlowcell(config, parkour, jobtrack.JobRegistry())

        assert set(secondPass) == set(dirs) - finished
        assert "beta" in secondPass
        for name in finished:
            assert name not in secondPass
        assert len(msg) == 3

    def test_partial_cleanup_leaves_a_corrupt_marker_undispatched(
        self, flowcell, monkeypatch
    ):
        config, parkour, dirs, _items = flowcell
        monkeypatch.setattr(jobtrack.subprocess, "run", lambda *a, **k: None)
        monkeypatch.setattr(PushButton, "POOL_SIZE", 1)

        firstPass = []
        registry1 = jobtrack.JobRegistry()
        gammaStarted = threading.Event()
        _stubPipeline(
            dirs,
            firstPass,
            {"beta": "crash", "gamma": "wait-for-abort"},
            monkeypatch,
            registry=registry1,
            started={"gamma": gammaStarted},
        )
        with pytest.raises(RuntimeError):
            PushButton.runFlowcell(config, parkour, registry1)
        # gamma is the "wait-for-abort" group -- settle its fate before our
        # manual corruption below, or it could race a still-running
        # straggler thread's own marker cleanup and get clobbered.
        _settleCrashedGroup(dirs, "gamma", gammaStarted)

        # The crash handler died between markMarkerCancelled and clearMarker,
        # leaving a truncated marker behind.
        (dirs["gamma"] / "running.pid").write_text('{"pid": 12, "job_i')

        secondPass = []
        _stubPipeline(dirs, secondPass, {}, monkeypatch)
        PushButton.runFlowcell(config, parkour, jobtrack.JobRegistry())

        assert "gamma" not in secondPass
        assert (dirs["gamma"] / "running.pid").exists()
        assert "beta" in secondPass

    def test_readable_cancelled_marker_from_a_dead_owner_is_redone(
        self, flowcell, monkeypatch
    ):
        config, parkour, dirs, _items = flowcell
        monkeypatch.setattr(jobtrack.subprocess, "run", lambda *a, **k: None)
        monkeypatch.setattr(PushButton, "POOL_SIZE", 1)
        jobtrack.writeMarker(dirs["gamma"], pid=999999, job_ids=["901"], cancelled=True)
        monkeypatch.setattr(jobtrack.os, "kill", _deadPid(999999))

        dispatched = []
        _stubPipeline(dirs, dispatched, {}, monkeypatch)
        PushButton.runFlowcell(config, parkour, jobtrack.JobRegistry())

        assert "gamma" in dispatched
        assert not (dirs["gamma"] / "running.pid").exists()

    def test_a_live_owner_blocks_redispatch(self, flowcell, monkeypatch):
        config, parkour, dirs, _items = flowcell
        monkeypatch.setattr(jobtrack.subprocess, "run", lambda *a, **k: None)
        monkeypatch.setattr(PushButton, "POOL_SIZE", 1)
        jobtrack.writeMarker(dirs["gamma"], job_ids=["901"])  # our own live pid

        dispatched = []
        _stubPipeline(dirs, dispatched, {}, monkeypatch)
        msg = PushButton.runFlowcell(config, parkour, jobtrack.JobRegistry())

        assert "gamma" not in dispatched
        assert set(dispatched) == {"alpha", "beta"}
        # Per Ruling 1 (Phase 1, commit 9e28c77), a skipped group still emits
        # a visible SKIPPED row rather than being silently dropped -- so all
        # three groups (2 dispatched + 1 skipped) show up in msg.
        assert len(msg) == 3
        statuses = {row[0]: row[4] for row in msg}
        assert statuses["2_Roe_Jones"] == "SKIPPED (owned by live PID)"


def _deadPid(dead):
    def _kill(pid, sig):
        if pid == dead:
            raise ProcessLookupError

    return _kill
