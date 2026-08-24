import configparser
import os
import threading
import time
from typing import ClassVar

import pytest

import BRB.ET
from BRB import PushButton, jobtrack
from BRB.PushButton import WorkItem

HUMAN = ("human", "hg38", "/yaml/hg38.yaml")
MOUSE = ("mouse", "mm10", "/yaml/mm10.yaml")


def dispatchConfig(tmp_path):
    """
    Config for the dispatch tests. groupData/<group> deliberately does NOT
    exist, so GetResults sets ignore=True and the [external] allowlist keeps
    the libraries in play -- which also makes createPath() root its output
    under baseData, with no PI group directory needed.
    """
    config = configparser.ConfigParser()
    config["Paths"] = {
        "baseData": str(tmp_path / "base"),
        "groupData": str(tmp_path / "group"),
        "seqFacDir": str(tmp_path / "seqfac"),
        "bioinfoCoreDir": str(tmp_path / "bioinfocore"),
    }
    config["Options"] = {
        "runID": "run1",
        "sequencerType": "Aviti",
        "validLibraryTypes": "ChIP-Seq,stranded mRNA-Seq",
        "pipelines": "DNA,RNA",
        "snakemakeWorkflowBaseDir": "/snakepipes",
    }
    config["external"] = {
        "LibraryTypes": "ChIP-Seq,stranded mRNA-Seq",
        "LibraryProtocols": "",
    }
    (tmp_path / "base" / "run1").mkdir(parents=True)
    (tmp_path / "group").mkdir(parents=True)
    (tmp_path / "bioinfocore").mkdir(parents=True)
    return config


class TestGetResultsReturnsWorkItems:
    def test_returns_workitems_and_does_not_dispatch(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)

        def boom(*args, **kwargs):
            raise AssertionError("GetResults must not dispatch pipelines any more")

        monkeypatch.setattr(PushButton, "DNA", boom)
        monkeypatch.setattr(PushButton, "RNA", boom)
        libraries = {
            "L1": ["s1", "ChIP-Seq", "proto", HUMAN, "i7", 30],
            "L2": ["s2", "stranded mRNA-Seq", "proto", HUMAN, "i7", 30],
        }

        workItems, msg = PushButton.GetResults(config, "1_A_Foo", libraries)

        assert msg == []
        assert len(workItems) == 2
        assert {(w.pipeline, w.libraryType) for w in workItems} == {
            ("DNA", "ChIP-Seq"),
            ("RNA", "stranded mRNA-Seq"),
        }
        for w in workItems:
            assert isinstance(w, WorkItem)
            assert w.project == "1_A_Foo"
            assert w.group == "foo"
            assert w.organism == HUMAN
        dnaItem = next(w for w in workItems if w.pipeline == "DNA")
        assert dnaItem.tuples == [["L1", "s1", "proto", True]]

    def test_skiplist_message_still_produced(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)
        # group dir exists => ignore=False => the invalid library type lands
        # in skipList rather than external_skipList.
        (tmp_path / "group" / "foo" / "sequencing_data").mkdir(parents=True)

        def fakeTelegraph(config, group, project, skipList, organism=None):
            return ["telegraph", organism, skipList[0][2], None, None, None, False]

        monkeypatch.setattr(BRB.ET, "telegraphHome", fakeTelegraph)
        libraries = {
            "L1": ["s1", "ChIP-Seq", "proto", HUMAN, "i7", 30],
            "L2": ["s2", "Other", "proto", HUMAN, "i7", 30],
        }

        workItems, msg = PushButton.GetResults(config, "1_A_Foo", libraries)

        assert len(workItems) == 1
        assert workItems[0].pipeline == "DNA"
        assert workItems[0].tuples == [["L1", "s1", "proto", False]]
        assert msg == [["telegraph", "human", "Other", None, None, None, False]]

    def test_external_skiplist_message_still_produced(self, tmp_path):
        config = dispatchConfig(tmp_path)
        config["external"]["LibraryTypes"] = "ChIP-Seq"
        libraries = {
            "L1": ["s1", "ChIP-Seq", "proto", HUMAN, "i7", 30],
            "L2": ["s2", "stranded mRNA-Seq", "proto", HUMAN, "i7", 30],
        }

        workItems, msg = PushButton.GetResults(config, "1_A_Foo", libraries)

        assert [w.libraryType for w in workItems] == ["ChIP-Seq"]
        assert msg == [
            ["1_A_Foo", "human", "stranded mRNA-Seq", None, None, None, False, None]
        ]


def makeWorkItem(**kwargs):
    base = {
        "project": "1_A_Foo",
        "group": "foo",
        "pipeline": "RNA",
        "organism": HUMAN,
        "libraryType": "stranded mRNA-Seq",
        "tuples": [["L1", "s1", "proto", True]],
    }
    base.update(kwargs)
    return WorkItem(**base)


def fakePhoneHome(config, outputDir, pipeline, tuples, organism, project, libType):
    return [project, organism, libType, pipeline, "success", "PARKOUR_OK"]


class TestRunOneGroup:
    def test_success_first_try(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)
        calls = []

        def stub(config, group, project, organism, libraryType, tuples):
            calls.append((group, project, libraryType))
            return "/out", 0, False

        monkeypatch.setattr(PushButton, "RNA", stub)
        monkeypatch.setattr(BRB.ET, "phoneHome", fakePhoneHome)

        result = PushButton.runOneGroup(config, makeWorkItem(), jobtrack.JobRegistry())

        assert len(calls) == 1
        assert calls[0] == ("foo", "1_A_Foo", "stranded mRNA-Seq")
        assert result == [
            [
                "1_A_Foo",
                "human",
                "stranded mRNA-Seq",
                "RNA",
                "success",
                "PARKOUR_OK",
                False,
                0,
            ]
        ]

    def test_retries_once_then_succeeds(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)
        calls = []

        def stub(config, group, project, organism, libraryType, tuples):
            calls.append(1)
            return ("/out", 1, False) if len(calls) == 1 else ("/out", 0, True)

        monkeypatch.setattr(PushButton, "RNA", stub)
        monkeypatch.setattr(BRB.ET, "phoneHome", fakePhoneHome)

        result = PushButton.runOneGroup(config, makeWorkItem(), jobtrack.JobRegistry())

        assert len(calls) == 2
        assert result[0][-2:] == [True, 1]
        assert result[0][4] == "success"

    def test_fails_twice_gives_FAILED_entry(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)
        calls = []

        def stub(config, group, project, organism, libraryType, tuples):
            calls.append(1)
            return "/out", 1, False

        monkeypatch.setattr(PushButton, "RNA", stub)

        def mustNotBeCalled(*args, **kwargs):
            raise AssertionError("phoneHome must not run for a failed group")

        monkeypatch.setattr(BRB.ET, "phoneHome", mustNotBeCalled)

        result = PushButton.runOneGroup(config, makeWorkItem(), jobtrack.JobRegistry())

        assert len(calls) == 2
        assert result == [
            [
                "1_A_Foo",
                "human",
                "stranded mRNA-Seq",
                "RNA",
                "FAILED",
                "not updated",
                False,
                1,
            ]
        ]


def findDeadPid():
    """A PID that is guaranteed not to name a live process right now."""
    pid = 4000000
    while pid > 2:
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return pid
        except (PermissionError, OverflowError, ValueError):
            pass
        pid -= 1
    raise RuntimeError("could not find a dead PID")


class TestOwnershipMarker:
    def _outputDir(self, config, item):
        return PushButton.createPath(
            config,
            item.group,
            item.project,
            item.organism[1],
            item.libraryType,
            item.tuples,
        )

    def test_no_marker_dispatches_and_cleans_up(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)
        item = makeWorkItem()
        calls = []
        outputDir = self._outputDir(config, item)

        def stub(config, group, project, organism, libraryType, tuples):
            calls.append(1)
            return outputDir, 0, False

        monkeypatch.setattr(PushButton, "RNA", stub)
        monkeypatch.setattr(BRB.ET, "phoneHome", fakePhoneHome)

        result = PushButton.runOneGroup(config, item, jobtrack.JobRegistry())

        assert len(calls) == 1
        assert len(result) == 1
        marker = os.path.join(outputDir, PushButton.MARKER_NAME)
        assert not os.path.exists(marker)

    def test_marker_written_while_dispatching(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)
        item = makeWorkItem()
        outputDir = self._outputDir(config, item)
        seen = {}

        def stub(config, group, project, organism, libraryType, tuples):
            seen["state"], marker = jobtrack.markerState(outputDir)
            seen["pid"] = marker["pid"]
            return outputDir, 0, False

        monkeypatch.setattr(PushButton, "RNA", stub)
        monkeypatch.setattr(BRB.ET, "phoneHome", fakePhoneHome)

        PushButton.runOneGroup(config, item, jobtrack.JobRegistry())

        assert seen["state"] == jobtrack.MARKER_LIVE
        assert seen["pid"] == os.getpid()

    def test_live_pid_marker_skips_group(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)
        item = makeWorkItem()
        outputDir = self._outputDir(config, item)
        marker = os.path.join(outputDir, PushButton.MARKER_NAME)
        jobtrack.writeMarker(outputDir)

        def mustNotBeCalled(*args, **kwargs):
            raise AssertionError("a live-PID-owned group must not be dispatched")

        monkeypatch.setattr(PushButton, "RNA", mustNotBeCalled)

        result = PushButton.runOneGroup(config, item, jobtrack.JobRegistry())

        assert result == [
            [
                "1_A_Foo",
                "human",
                "stranded mRNA-Seq",
                "RNA",
                "SKIPPED (owned by live PID)",
                "not updated",
                False,
                0,
            ]
        ]
        assert os.path.exists(marker), "another owner's marker must not be removed"

    def test_dead_pid_marker_is_abandoned_and_dispatch_proceeds(
        self, tmp_path, monkeypatch
    ):
        config = dispatchConfig(tmp_path)
        item = makeWorkItem()
        outputDir = self._outputDir(config, item)
        marker = os.path.join(outputDir, PushButton.MARKER_NAME)
        jobtrack.writeMarker(outputDir, pid=findDeadPid())
        calls = []

        def stub(config, group, project, organism, libraryType, tuples):
            calls.append(1)
            return outputDir, 0, False

        monkeypatch.setattr(PushButton, "RNA", stub)
        monkeypatch.setattr(BRB.ET, "phoneHome", fakePhoneHome)

        result = PushButton.runOneGroup(config, item, jobtrack.JobRegistry())

        assert len(calls) == 1
        assert len(result) == 1
        assert not os.path.exists(marker)

    def test_corrupt_marker_skips_group_and_leaves_marker_in_place(
        self, tmp_path, monkeypatch
    ):
        """
        Superseded/replaces the old non-positive-pid case: jobtrack's marker
        is a JSON dict, so a marker written with a non-integer/missing "pid"
        field is CORRUPT (unreadable), not a special-cased "not live" PID --
        it is skipped, not dispatched, exactly like a truly unreadable file.
        """
        config = dispatchConfig(tmp_path)
        item = makeWorkItem()
        outputDir = self._outputDir(config, item)
        marker = os.path.join(outputDir, PushButton.MARKER_NAME)
        with open(marker, "w") as fh:
            fh.write('{"pid": "not-a-pid"}')

        def mustNotBeCalled(*args, **kwargs):
            raise AssertionError("a corrupt-marker-owned group must not be dispatched")

        monkeypatch.setattr(PushButton, "RNA", mustNotBeCalled)

        result = PushButton.runOneGroup(config, item, jobtrack.JobRegistry())

        assert result[0][4] == "SKIPPED (unreadable marker)"
        assert os.path.exists(marker), "an unreadable marker must not be removed"

    def test_marker_removed_when_pipeline_raises(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)
        item = makeWorkItem()

        def stub(config, group, project, organism, libraryType, tuples):
            raise RuntimeError("boom")

        monkeypatch.setattr(PushButton, "RNA", stub)

        with pytest.raises(RuntimeError):
            PushButton.runOneGroup(config, item, jobtrack.JobRegistry())

        marker = os.path.join(self._outputDir(config, item), PushButton.MARKER_NAME)
        assert not os.path.exists(marker)


class TestRunFlowcell:
    def _parkourDict(self):
        return {
            "1_A_Foo": {
                "L1": ["s1", "ChIP-Seq", "proto", HUMAN, "i7", 30],
                "L2": ["s2", "stranded mRNA-Seq", "proto", HUMAN, "i7", 30],
            },
            "2_B_Bar": {"L3": ["s3", "ChIP-Seq", "proto", MOUSE, "i7", 30]},
        }

    def _makeProjectDirs(self, tmp_path, names):
        for name in names:
            (tmp_path / "base" / "run1" / f"Project_{name}").mkdir(parents=True)

    def test_collects_messages_from_every_group(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)
        self._makeProjectDirs(tmp_path, ["1_A_Foo", "2_B_Bar"])
        monkeypatch.setattr(PushButton, "POOL_SIZE", 2)
        monkeypatch.setattr(BRB.ET, "phoneHome", fakePhoneHome)

        def stub(config, group, project, organism, libraryType, tuples):
            return "/out", 0, False

        monkeypatch.setattr(PushButton, "DNA", stub)
        monkeypatch.setattr(PushButton, "RNA", stub)

        msg = PushButton.runFlowcell(config, self._parkourDict())

        assert sorted((m[0], m[2], m[3]) for m in msg) == [
            ("1_A_Foo", "ChIP-Seq", "DNA"),
            ("1_A_Foo", "stranded mRNA-Seq", "RNA"),
            ("2_B_Bar", "ChIP-Seq", "DNA"),
        ]

    def test_skips_project_that_lives_on_another_lane(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)
        self._makeProjectDirs(tmp_path, ["2_B_Bar"])  # 1_A_Foo absent
        monkeypatch.setattr(PushButton, "POOL_SIZE", 2)
        monkeypatch.setattr(BRB.ET, "phoneHome", fakePhoneHome)

        def stub(config, group, project, organism, libraryType, tuples):
            return "/out", 0, False

        monkeypatch.setattr(PushButton, "DNA", stub)
        monkeypatch.setattr(PushButton, "RNA", stub)

        msg = PushButton.runFlowcell(config, self._parkourDict())

        assert [m[0] for m in msg] == ["2_B_Bar"]

    def test_no_work_items_returns_empty(self, tmp_path, monkeypatch):
        config = dispatchConfig(tmp_path)
        monkeypatch.setattr(PushButton, "POOL_SIZE", 2)
        assert PushButton.runFlowcell(config, {}) == []

    def test_marker_owned_group_appears_as_skipped_not_silently_absent(
        self, tmp_path, monkeypatch
    ):
        """
        runFlowcell-level (not just runOneGroup-level) proof that a
        marker-owned group shows up as a visible SKIPPED entry in the
        returned msg, alongside the other groups' normal results -- it must
        not be silently absent, which would otherwise let run.py's
        subsequent markFinished() call mark the flowcell fully done even
        though this group was never analysed this pass.
        """
        config = dispatchConfig(tmp_path)
        self._makeProjectDirs(tmp_path, ["1_A_Foo", "2_B_Bar"])
        monkeypatch.setattr(PushButton, "POOL_SIZE", 2)
        monkeypatch.setattr(BRB.ET, "phoneHome", fakePhoneHome)

        # Pre-seed a live-PID ownership marker for 1_A_Foo's ChIP-Seq (DNA)
        # group, so runOneGroup skips it instead of dispatching.
        outputDir = PushButton.createPath(
            config, "foo", "1_A_Foo", "hg38", "ChIP-Seq", [[None, None, None, True]]
        )
        marker = os.path.join(outputDir, PushButton.MARKER_NAME)
        jobtrack.writeMarker(outputDir)

        def dnaStub(config, group, project, organism, libraryType, tuples):
            assert not (project == "1_A_Foo" and libraryType == "ChIP-Seq"), (
                "the marker-owned group must not be dispatched"
            )
            return "/out", 0, False

        def rnaStub(config, group, project, organism, libraryType, tuples):
            return "/out", 0, False

        monkeypatch.setattr(PushButton, "DNA", dnaStub)
        monkeypatch.setattr(PushButton, "RNA", rnaStub)

        msg = PushButton.runFlowcell(config, self._parkourDict())

        skipped = [m for m in msg if m[0] == "1_A_Foo" and m[2] == "ChIP-Seq"]
        assert len(skipped) == 1, f"marker-owned group missing from msg: {msg}"
        assert skipped[0][4] == "SKIPPED (owned by live PID)"
        others = [m for m in msg if not (m[0] == "1_A_Foo" and m[2] == "ChIP-Seq")]
        assert sorted((m[0], m[2]) for m in others) == [
            ("1_A_Foo", "stranded mRNA-Seq"),
            ("2_B_Bar", "ChIP-Seq"),
        ]
        assert os.path.exists(marker), "another owner's marker must not be removed"

    def test_crash_path_does_not_block_on_surviving_workers(
        self, tmp_path, monkeypatch
    ):
        config = dispatchConfig(tmp_path)
        self._makeProjectDirs(tmp_path, ["1_A_Foo", "2_B_Bar"])
        monkeypatch.setattr(PushButton, "POOL_SIZE", 2)
        monkeypatch.setattr(BRB.ET, "phoneHome", fakePhoneHome)

        def slowStub(config, group, project, organism, libraryType, tuples):
            time.sleep(3)
            return "/out", 0, False

        def crashStub(config, group, project, organism, libraryType, tuples):
            raise RuntimeError("kaboom")

        # DNA groups sleep, the RNA group crashes immediately.
        monkeypatch.setattr(PushButton, "DNA", slowStub)
        monkeypatch.setattr(PushButton, "RNA", crashStub)

        started = time.monotonic()
        with pytest.raises(PushButton.GroupDispatchError) as excinfo:
            PushButton.runFlowcell(config, self._parkourDict())
        elapsed = time.monotonic() - started

        assert elapsed < 1.5, "crash path must not wait for surviving workers"
        assert isinstance(excinfo.value.__cause__, RuntimeError)
        failedItems = [item for item, _exc in excinfo.value.failures]
        assert any(
            i.pipeline == "RNA" and i.libraryType == "stranded mRNA-Seq"
            for i in failedItems
        )
        assert "1_A_Foo" in str(excinfo.value)
        assert "RNA" in str(excinfo.value)
        assert "stranded mRNA-Seq" in str(excinfo.value)

    def test_shutdown_called_with_wait_false_and_cancel_futures(
        self, tmp_path, monkeypatch
    ):
        config = dispatchConfig(tmp_path)
        self._makeProjectDirs(tmp_path, ["1_A_Foo", "2_B_Bar"])
        monkeypatch.setattr(PushButton, "POOL_SIZE", 2)
        recorded = {}

        realExecutor = PushButton.ThreadPoolExecutor

        class RecordingExecutor(realExecutor):
            def shutdown(self, wait=True, *, cancel_futures=False):
                recorded.setdefault("calls", []).append((wait, cancel_futures))
                return super().shutdown(wait=False, cancel_futures=cancel_futures)

        monkeypatch.setattr(PushButton, "ThreadPoolExecutor", RecordingExecutor)

        def crashStub(config, group, project, organism, libraryType, tuples):
            raise RuntimeError("kaboom")

        monkeypatch.setattr(PushButton, "DNA", crashStub)
        monkeypatch.setattr(PushButton, "RNA", crashStub)

        with pytest.raises(PushButton.GroupDispatchError):
            PushButton.runFlowcell(config, self._parkourDict())

        assert (False, True) in recorded["calls"]

    def test_every_collected_failure_is_logged_critical(
        self, tmp_path, monkeypatch, caplog
    ):
        """
        The invariant under test is that runFlowcell logs one log.critical
        line per entry in the raised GroupDispatchError.failures, naming that
        entry's project/pipeline/libraryType -- i.e. it never reports only the
        first failure it saw.

        Deliberately NOT asserting `len(failures) >= 2`: how many sibling
        futures have already reached `done()` at the instant the first
        exception surfaces is genuinely scheduler-dependent, and pinning it
        would make this test flaky rather than stricter.
        """
        config = dispatchConfig(tmp_path)
        self._makeProjectDirs(tmp_path, ["1_A_Foo", "2_B_Bar"])
        monkeypatch.setattr(PushButton, "POOL_SIZE", 4)

        def crashStub(config, group, project, organism, libraryType, tuples):
            raise RuntimeError(f"kaboom {project} {libraryType}")

        monkeypatch.setattr(PushButton, "DNA", crashStub)
        monkeypatch.setattr(PushButton, "RNA", crashStub)

        with (
            caplog.at_level("CRITICAL"),
            pytest.raises(PushButton.GroupDispatchError) as excinfo,
        ):
            PushButton.runFlowcell(config, self._parkourDict())

        criticals = [
            r.getMessage() for r in caplog.records if r.levelname == "CRITICAL"
        ]
        failures = excinfo.value.failures
        assert len(failures) >= 1
        assert len(criticals) == len(failures)
        for item, _exc in failures:
            assert any(
                item.project in m and item.pipeline in m and item.libraryType in m
                for m in criticals
            ), f"no log.critical line for {item.project}/{item.pipeline}"
            assert item.project in str(excinfo.value)


class ConcurrencyTracker:
    """Shared across the stub pipelines so the cap is measured pool-wide."""

    def __init__(self):
        self.lock = threading.Lock()
        self.active = 0
        self.maxActive = 0
        self.calls = []

    def enter(self, project, libraryType):
        with self.lock:
            self.active += 1
            self.maxActive = max(self.maxActive, self.active)
            self.calls.append((project, libraryType))

    def leave(self):
        with self.lock:
            self.active -= 1


def makeStub(tracker, rv, delay=0.2, sambaUpdate=False):
    def stub(config, group, project, organism, libraryType, tuples):
        tracker.enter(project, libraryType)
        try:
            time.sleep(delay)
        finally:
            tracker.leave()
        # Return the *real* outputDir (matching what every actual pipeline
        # function does: it calls createPath with the same arguments it was
        # given, which is deterministic). runOneGroup's marker cleanup and
        # registry.unregister_group() key off the value returned here, so a
        # fake/mismatched path would leave the real marker file behind.
        outputDir = PushButton.createPath(
            config, group, project, organism[1], libraryType, tuples
        )
        return outputDir, rv, sambaUpdate

    return stub


class TestRunFlowcellIntegration:
    PARKOUR: ClassVar = {
        "1_A_Foo": {
            "L1": ["s1", "ChIP-Seq", "proto", HUMAN, "i7", 30],
            "L2": ["s2", "stranded mRNA-Seq", "proto", HUMAN, "i7", 30],
        },
        "2_B_Bar": {"L3": ["s3", "ChIP-Seq", "proto", MOUSE, "i7", 30]},
        "3_C_Baz": {"L4": ["s4", "stranded mRNA-Seq", "proto", MOUSE, "i7", 30]},
    }

    def _setup(self, tmp_path, monkeypatch, poolSize=2):
        config = dispatchConfig(tmp_path)
        for name in self.PARKOUR:
            (tmp_path / "base" / "run1" / f"Project_{name}").mkdir(parents=True)
        monkeypatch.setattr(PushButton, "POOL_SIZE", poolSize)
        monkeypatch.setattr(BRB.ET, "phoneHome", fakePhoneHome)
        return config

    def test_pool_size_is_respected_and_all_groups_dispatched(
        self, tmp_path, monkeypatch
    ):
        config = self._setup(tmp_path, monkeypatch, poolSize=2)
        tracker = ConcurrencyTracker()
        monkeypatch.setattr(PushButton, "DNA", makeStub(tracker, rv=0))
        monkeypatch.setattr(PushButton, "RNA", makeStub(tracker, rv=0))

        msg = PushButton.runFlowcell(config, self.PARKOUR)

        assert tracker.maxActive <= 2, "more than POOL_SIZE groups ran at once"
        assert tracker.maxActive == 2, "the pool never actually overlapped groups"
        assert sorted(tracker.calls) == sorted(
            [
                ("1_A_Foo", "ChIP-Seq"),
                ("1_A_Foo", "stranded mRNA-Seq"),
                ("2_B_Bar", "ChIP-Seq"),
                ("3_C_Baz", "stranded mRNA-Seq"),
            ]
        )
        assert len(msg) == 4

    def test_results_collect_order_independently_and_failures_stay_local(
        self, tmp_path, monkeypatch
    ):
        config = self._setup(tmp_path, monkeypatch, poolSize=2)
        tracker = ConcurrencyTracker()
        # DNA groups succeed quickly; RNA groups always fail, so each is
        # retried once and yields a FAILED entry -- without blocking or
        # failing the DNA groups.
        monkeypatch.setattr(PushButton, "DNA", makeStub(tracker, rv=0, delay=0.3))
        monkeypatch.setattr(PushButton, "RNA", makeStub(tracker, rv=1, delay=0.05))

        msg = PushButton.runFlowcell(config, self.PARKOUR)

        assert sorted(tuple(m) for m in msg) == sorted(
            [
                (
                    "1_A_Foo",
                    "human",
                    "ChIP-Seq",
                    "DNA",
                    "success",
                    "PARKOUR_OK",
                    False,
                    0,
                ),
                (
                    "2_B_Bar",
                    "mouse",
                    "ChIP-Seq",
                    "DNA",
                    "success",
                    "PARKOUR_OK",
                    False,
                    0,
                ),
                (
                    "1_A_Foo",
                    "human",
                    "stranded mRNA-Seq",
                    "RNA",
                    "FAILED",
                    "not updated",
                    False,
                    1,
                ),
                (
                    "3_C_Baz",
                    "mouse",
                    "stranded mRNA-Seq",
                    "RNA",
                    "FAILED",
                    "not updated",
                    False,
                    1,
                ),
            ]
        )
        # Each RNA group was attempted twice (initial + one retry).
        rnaCalls = [c for c in tracker.calls if c[1] == "stranded mRNA-Seq"]
        assert len(rnaCalls) == 4

    def test_no_ownership_markers_left_behind(self, tmp_path, monkeypatch):
        config = self._setup(tmp_path, monkeypatch, poolSize=2)
        tracker = ConcurrencyTracker()
        monkeypatch.setattr(PushButton, "DNA", makeStub(tracker, rv=0, delay=0.05))
        monkeypatch.setattr(PushButton, "RNA", makeStub(tracker, rv=0, delay=0.05))

        PushButton.runFlowcell(config, self.PARKOUR)

        leftovers = list((tmp_path / "base").rglob(PushButton.MARKER_NAME))
        assert leftovers == []
