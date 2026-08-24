import json
import os
import threading
import time

import pytest

from BRB import jobtrack


class TestJobRegistry:
    def test_register_returns_handle_for_outputdir(self, tmp_path):
        reg = jobtrack.JobRegistry()
        handle = reg.register_group(tmp_path)
        assert handle.outputDir == str(tmp_path)

    def test_register_is_idempotent_per_outputdir(self, tmp_path):
        reg = jobtrack.JobRegistry()
        assert reg.register_group(tmp_path) is reg.register_group(str(tmp_path))

    def test_active_groups_lists_registered_handles(self, tmp_path):
        reg = jobtrack.JobRegistry()
        a = reg.register_group(tmp_path / "a")
        b = reg.register_group(tmp_path / "b")
        assert set(reg.active_groups()) == {a, b}

    def test_unregister_removes_from_active(self, tmp_path):
        reg = jobtrack.JobRegistry()
        handle = reg.register_group(tmp_path)
        assert reg.unregister_group(tmp_path) is handle
        assert reg.active_groups() == []

    def test_unregister_unknown_is_none(self, tmp_path):
        assert jobtrack.JobRegistry().unregister_group(tmp_path) is None

    def test_two_registries_are_independent(self, tmp_path):
        one, two = jobtrack.JobRegistry(), jobtrack.JobRegistry()
        one.register_group(tmp_path)
        assert two.active_groups() == []

    def test_abort_flag_defaults_false_and_latches(self):
        reg = jobtrack.JobRegistry()
        assert reg.aborted is False
        reg.abort()
        assert reg.aborted is True

    def test_concurrent_registration_is_consistent(self, tmp_path):
        reg = jobtrack.JobRegistry()
        seen = []
        barrier = threading.Barrier(8)

        def worker():
            barrier.wait()
            seen.append(reg.register_group(tmp_path))

        threads = [threading.Thread(target=worker) for _ in range(8)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        assert len(set(seen)) == 1
        assert len(reg.active_groups()) == 1


class TestGroupHandle:
    def test_new_handle_has_no_jobs_or_processes(self, tmp_path):
        handle = jobtrack.JobRegistry().register_group(tmp_path)
        assert handle.snapshot() == ([], [])

    def test_add_job_ids_dedupes_and_preserves_order(self, tmp_path):
        handle = jobtrack.JobRegistry().register_group(tmp_path)
        assert handle.add_job_ids(["101", "102"]) == ["101", "102"]
        assert handle.add_job_ids(["102", "103"]) == ["103"]
        assert handle.snapshot()[0] == ["101", "102", "103"]

    def test_process_add_and_remove(self, tmp_path):
        handle = jobtrack.JobRegistry().register_group(tmp_path)
        proc = object()
        handle.add_process(proc)
        assert handle.snapshot()[1] == [proc]
        handle.remove_process(proc)
        assert handle.snapshot()[1] == []

    def test_remove_unknown_process_is_a_noop(self, tmp_path):
        handle = jobtrack.JobRegistry().register_group(tmp_path)
        handle.remove_process(object())
        assert handle.snapshot()[1] == []

    def test_snapshot_returns_copies(self, tmp_path):
        handle = jobtrack.JobRegistry().register_group(tmp_path)
        handle.add_job_ids(["1"])
        jobs, _ = handle.snapshot()
        jobs.append("2")
        assert handle.snapshot()[0] == ["1"]


class TestBindGroup:
    def test_current_handle_is_none_outside_a_binding(self):
        assert jobtrack.currentHandle() is None

    def test_bind_group_sets_and_restores(self, tmp_path):
        handle = jobtrack.JobRegistry().register_group(tmp_path)
        with jobtrack.bindGroup(handle) as bound:
            assert bound is handle
            assert jobtrack.currentHandle() is handle
        assert jobtrack.currentHandle() is None

    def test_bind_group_restores_on_exception(self, tmp_path):
        handle = jobtrack.JobRegistry().register_group(tmp_path)
        with pytest.raises(RuntimeError), jobtrack.bindGroup(handle):
            raise RuntimeError("boom")
        assert jobtrack.currentHandle() is None

    def test_binding_is_per_thread(self, tmp_path):
        reg = jobtrack.JobRegistry()
        outer = reg.register_group(tmp_path / "outer")
        seen = []

        def worker():
            seen.append(jobtrack.currentHandle())

        with jobtrack.bindGroup(outer):
            t = threading.Thread(target=worker)
            t.start()
            t.join()
        assert seen == [None]


class TestMarker:
    def test_write_then_read_round_trips(self, tmp_path):
        written = jobtrack.writeMarker(tmp_path, job_ids=["7", "8"])
        read = jobtrack.readMarker(tmp_path)
        assert read == written
        assert read["pid"] == os.getpid()
        assert read["job_ids"] == ["7", "8"]
        assert read["cancelled"] is False
        assert isinstance(read["started"], float)

    def test_marker_lands_at_running_pid(self, tmp_path):
        jobtrack.writeMarker(tmp_path)
        assert (tmp_path / "running.pid").exists()
        assert jobtrack.MARKER_NAME == "running.pid"

    def test_write_defaults_to_empty_job_ids(self, tmp_path):
        assert jobtrack.writeMarker(tmp_path)["job_ids"] == []

    def test_read_absent_marker_is_none(self, tmp_path):
        assert jobtrack.readMarker(tmp_path) is None

    def test_update_job_ids_preserves_pid_and_started(self, tmp_path):
        first = jobtrack.writeMarker(tmp_path, pid=4242)
        updated = jobtrack.updateMarkerJobIds(tmp_path, ["101"])
        assert updated["pid"] == 4242
        assert updated["started"] == first["started"]
        assert updated["job_ids"] == ["101"]

    def test_update_job_ids_on_absent_marker_creates_one(self, tmp_path):
        updated = jobtrack.updateMarkerJobIds(tmp_path, ["55"])
        assert jobtrack.readMarker(tmp_path)["job_ids"] == ["55"]
        assert updated["pid"] == os.getpid()

    def test_mark_cancelled_sets_the_flag(self, tmp_path):
        jobtrack.writeMarker(tmp_path, job_ids=["9"])
        assert jobtrack.markMarkerCancelled(tmp_path)["cancelled"] is True
        assert jobtrack.readMarker(tmp_path)["job_ids"] == ["9"]

    def test_mark_cancelled_on_absent_marker_is_none(self, tmp_path):
        assert jobtrack.markMarkerCancelled(tmp_path) is None

    def test_clear_marker_removes_it(self, tmp_path):
        jobtrack.writeMarker(tmp_path)
        jobtrack.clearMarker(tmp_path)
        assert not (tmp_path / "running.pid").exists()

    def test_clear_marker_is_idempotent(self, tmp_path):
        jobtrack.clearMarker(tmp_path)
        jobtrack.clearMarker(tmp_path)

    def test_no_tempfile_left_behind(self, tmp_path):
        jobtrack.writeMarker(tmp_path, job_ids=["1"])
        jobtrack.updateMarkerJobIds(tmp_path, ["1", "2"])
        assert [p.name for p in tmp_path.iterdir()] == ["running.pid"]


class TestMarkerState:
    def test_absent(self, tmp_path):
        assert jobtrack.markerState(tmp_path) == (jobtrack.MARKER_ABSENT, None)

    def test_live_pid_is_live(self, tmp_path):
        jobtrack.writeMarker(tmp_path)
        state, marker = jobtrack.markerState(tmp_path)
        assert state == jobtrack.MARKER_LIVE
        assert marker["pid"] == os.getpid()

    def test_dead_pid_is_abandoned(self, tmp_path, monkeypatch):
        jobtrack.writeMarker(tmp_path, pid=999999)
        monkeypatch.setattr(jobtrack.os, "kill", _raise_no_such_process, raising=True)
        state, marker = jobtrack.markerState(tmp_path)
        assert state == jobtrack.MARKER_ABANDONED
        assert marker["pid"] == 999999

    def test_cancelled_dead_marker_is_still_abandoned(self, tmp_path, monkeypatch):
        jobtrack.writeMarker(tmp_path, pid=999999, job_ids=["3"], cancelled=True)
        monkeypatch.setattr(jobtrack.os, "kill", _raise_no_such_process)
        state, marker = jobtrack.markerState(tmp_path)
        assert state == jobtrack.MARKER_ABANDONED
        assert marker["cancelled"] is True

    def test_truncated_marker_is_corrupt(self, tmp_path):
        (tmp_path / "running.pid").write_text('{"pid": 12')
        assert jobtrack.markerState(tmp_path) == (jobtrack.MARKER_CORRUPT, None)

    def test_marker_without_pid_is_corrupt(self, tmp_path):
        (tmp_path / "running.pid").write_text(json.dumps({"job_ids": []}))
        assert jobtrack.markerState(tmp_path) == (jobtrack.MARKER_CORRUPT, None)

    def test_non_object_marker_is_corrupt(self, tmp_path):
        (tmp_path / "running.pid").write_text("[1, 2, 3]")
        assert jobtrack.markerState(tmp_path) == (jobtrack.MARKER_CORRUPT, None)

    def test_permission_error_counts_as_live(self, tmp_path, monkeypatch):
        jobtrack.writeMarker(tmp_path, pid=1)

        def _perm(pid, sig):
            raise PermissionError

        monkeypatch.setattr(jobtrack.os, "kill", _perm)
        assert jobtrack.markerState(tmp_path)[0] == jobtrack.MARKER_LIVE


def _raise_no_such_process(pid, sig):
    raise ProcessLookupError


class TestParseJobIds:
    def test_cluster_generic_submission_line(self):
        line = "Submitted job 12 with external jobid 'Submitted batch job 1234567'."
        assert jobtrack.parseJobIds(line) == ["1234567"]

    def test_parsable_sbatch_bare_id(self):
        assert jobtrack.parseJobIds("external jobid '7654321'") == ["7654321"]

    def test_job_id_with_cluster_suffix(self):
        line = "Submitted job 3 with external jobid 'Submitted batch job 42 on cluster deep'."
        assert jobtrack.parseJobIds(line) == ["42"]

    def test_unrelated_lines_yield_nothing(self):
        assert jobtrack.parseJobIds("rule bowtie2: output=foo.bam") == []
        assert jobtrack.parseJobIds("") == []
        assert jobtrack.parseJobIds("Job 12 has been submitted") == []

    def test_numbers_elsewhere_are_not_job_ids(self):
        assert jobtrack.parseJobIds("Provided cores: 16") == []

    def test_two_ids_on_one_line(self):
        line = "external jobid '111' and external jobid '222'"
        assert jobtrack.parseJobIds(line) == ["111", "222"]

    def test_empty_jobid_payload_is_ignored(self):
        assert jobtrack.parseJobIds("external jobid ''") == []

    def test_production_cluster_parsable_format(self):
        """MPI-IE cluster uses sbatch --parsable, producing bare integer IDs."""
        assert jobtrack.parseJobIds(
            "Submitted job 1 with external jobid '10874970'."
        ) == ["10874970"]


import subprocess


class TestRunManagedSubprocess:
    def test_returns_zero_on_success(self):
        assert jobtrack.runManagedSubprocess("true") == 0

    def test_raises_calledprocesserror_on_failure(self):
        with pytest.raises(subprocess.CalledProcessError) as exc:
            jobtrack.runManagedSubprocess("exit 3")
        assert exc.value.returncode == 3

    def test_runs_in_cwd_when_given(self, tmp_path):
        jobtrack.runManagedSubprocess("pwd > where.txt", cwd=tmp_path)
        assert (tmp_path / "where.txt").read_text().strip() == str(tmp_path)

    def test_child_gets_its_own_session(self):
        jobtrack.runManagedSubprocess('test "$(ps -o sid= -p $$ | tr -d \' \')" = "$$"')

    def test_output_is_logged(self, caplog):
        with caplog.at_level("INFO"):
            jobtrack.runManagedSubprocess("echo hello-from-driver")
        assert "hello-from-driver" in caplog.text

    def test_stderr_is_captured_too(self, caplog):
        with caplog.at_level("INFO"):
            jobtrack.runManagedSubprocess("echo oops >&2")
        assert "oops" in caplog.text

    def test_job_ids_land_in_the_bound_handle(self, tmp_path):
        handle = jobtrack.JobRegistry().register_group(tmp_path)
        with jobtrack.bindGroup(handle):
            jobtrack.runManagedSubprocess(
                'echo "Submitted job 1 with external jobid '
                "'Submitted batch job 555'.\""
            )
        assert handle.snapshot()[0] == ["555"]

    def test_job_ids_land_in_the_marker_file(self, tmp_path):
        handle = jobtrack.JobRegistry().register_group(tmp_path)
        jobtrack.writeMarker(tmp_path)
        with jobtrack.bindGroup(handle):
            jobtrack.runManagedSubprocess("echo \"external jobid '777'\"")
        assert jobtrack.readMarker(tmp_path)["job_ids"] == ["777"]

    def test_explicit_handle_beats_the_bound_one(self, tmp_path):
        reg = jobtrack.JobRegistry()
        bound = reg.register_group(tmp_path / "bound")
        explicit = reg.register_group(tmp_path / "explicit")
        (tmp_path / "explicit").mkdir()
        with jobtrack.bindGroup(bound):
            jobtrack.runManagedSubprocess(
                "echo \"external jobid '888'\"", handle=explicit
            )
        assert explicit.snapshot()[0] == ["888"]
        assert bound.snapshot()[0] == []

    def test_no_handle_bound_is_fine(self):
        assert jobtrack.runManagedSubprocess("echo \"external jobid '999'\"") == 0

    def test_process_is_deregistered_after_success(self, tmp_path):
        handle = jobtrack.JobRegistry().register_group(tmp_path)
        with jobtrack.bindGroup(handle):
            jobtrack.runManagedSubprocess("true")
        assert handle.snapshot()[1] == []

    def test_process_is_deregistered_after_failure(self, tmp_path):
        handle = jobtrack.JobRegistry().register_group(tmp_path)
        with jobtrack.bindGroup(handle), pytest.raises(subprocess.CalledProcessError):
            jobtrack.runManagedSubprocess("exit 1")
        assert handle.snapshot()[1] == []

    def test_process_is_visible_while_running(self, tmp_path):
        handle = jobtrack.JobRegistry().register_group(tmp_path)
        seen = []

        def watcher():
            for _ in range(200):
                if handle.snapshot()[1]:
                    seen.append(True)
                    return
                time.sleep(0.01)

        t = threading.Thread(target=watcher)
        t.start()
        with jobtrack.bindGroup(handle):
            jobtrack.runManagedSubprocess("sleep 0.5")
        t.join()
        assert seen == [True]


import signal


class FakeProc:
    def __init__(self, pid=4242):
        self.pid = pid
        self._exited = False

    def poll(self):
        return 0 if self._exited else None

    def wait(self, timeout=None):
        self._exited = True
        return 0


class TestCancelGroup:
    def test_scancels_tracked_job_ids(self, tmp_path, monkeypatch):
        handle = jobtrack.JobRegistry().register_group(tmp_path)
        handle.add_job_ids(["101", "102"])
        runs = []
        monkeypatch.setattr(
            jobtrack.subprocess, "run", lambda *a, **k: runs.append((a, k))
        )
        jobtrack.cancelGroup(handle)
        assert runs[0][0][0] == ["scancel", "101", "102"]

    def test_no_jobs_means_no_scancel(self, tmp_path, monkeypatch):
        handle = jobtrack.JobRegistry().register_group(tmp_path)
        runs = []
        monkeypatch.setattr(jobtrack.subprocess, "run", lambda *a, **k: runs.append(a))
        jobtrack.cancelGroup(handle)
        assert runs == []

    def test_sigterms_each_driver_process_group(self, tmp_path, monkeypatch):
        handle = jobtrack.JobRegistry().register_group(tmp_path)
        proc = FakeProc(pid=777)
        handle.add_process(proc)
        signals = []
        monkeypatch.setattr(jobtrack.os, "getpgid", lambda pid: pid)
        monkeypatch.setattr(
            jobtrack.os, "killpg", lambda pgid, sig: signals.append((pgid, sig))
        )
        jobtrack.cancelGroup(handle)
        assert signals == [(777, signal.SIGTERM)]

    def test_marker_is_flagged_cancelled_before_the_kill(self, tmp_path, monkeypatch):
        handle = jobtrack.JobRegistry().register_group(tmp_path)
        handle.add_job_ids(["55"])
        jobtrack.writeMarker(tmp_path, job_ids=["55"])
        order = []
        monkeypatch.setattr(
            jobtrack.subprocess,
            "run",
            lambda *a, **k: order.append(jobtrack.readMarker(tmp_path)["cancelled"]),
        )
        jobtrack.cancelGroup(handle)
        assert order == [True]

    def test_scancel_failure_does_not_stop_the_sigterm(self, tmp_path, monkeypatch):
        handle = jobtrack.JobRegistry().register_group(tmp_path)
        handle.add_job_ids(["1"])
        proc = FakeProc(pid=888)
        handle.add_process(proc)
        signals = []

        def blowUp(*a, **k):
            raise OSError("scancel: command not found")

        monkeypatch.setattr(jobtrack.subprocess, "run", blowUp)
        monkeypatch.setattr(jobtrack.os, "getpgid", lambda pid: pid)
        monkeypatch.setattr(
            jobtrack.os, "killpg", lambda pgid, sig: signals.append((pgid, sig))
        )
        jobtrack.cancelGroup(handle)
        assert signals == [(888, signal.SIGTERM)]

    def test_already_dead_process_is_tolerated(self, tmp_path, monkeypatch):
        handle = jobtrack.JobRegistry().register_group(tmp_path)
        handle.add_process(FakeProc(pid=999))

        def gone(pid):
            raise ProcessLookupError

        monkeypatch.setattr(jobtrack.os, "getpgid", gone)
        jobtrack.cancelGroup(handle)


class TestCancelAllGroups:
    def test_cancels_every_registered_group(self, tmp_path, monkeypatch):
        reg = jobtrack.JobRegistry()
        for name in ("a", "b", "c"):
            (tmp_path / name).mkdir()
            handle = reg.register_group(tmp_path / name)
            handle.add_job_ids([f"{name}1"])
            jobtrack.writeMarker(tmp_path / name)
        runs = []
        monkeypatch.setattr(
            jobtrack.subprocess, "run", lambda *a, **k: runs.append(a[0])
        )
        jobtrack.cancelAllGroups(reg)
        assert sorted(r[1] for r in runs) == ["a1", "b1", "c1"]

    def test_clears_markers_and_empties_the_registry(self, tmp_path, monkeypatch):
        reg = jobtrack.JobRegistry()
        (tmp_path / "a").mkdir()
        reg.register_group(tmp_path / "a")
        jobtrack.writeMarker(tmp_path / "a")
        monkeypatch.setattr(jobtrack.subprocess, "run", lambda *a, **k: None)
        jobtrack.cancelAllGroups(reg)
        assert not (tmp_path / "a" / "running.pid").exists()
        assert reg.active_groups() == []

    def test_sets_the_abort_flag(self, tmp_path, monkeypatch):
        reg = jobtrack.JobRegistry()
        monkeypatch.setattr(jobtrack.subprocess, "run", lambda *a, **k: None)
        jobtrack.cancelAllGroups(reg)
        assert reg.aborted is True

    def test_one_group_failing_does_not_stop_the_others(self, tmp_path, monkeypatch):
        reg = jobtrack.JobRegistry()
        for name in ("a", "b"):
            (tmp_path / name).mkdir()
            reg.register_group(tmp_path / name)
            jobtrack.writeMarker(tmp_path / name)
        calls = []

        def flaky(handle, **kw):
            calls.append(handle.outputDir)
            if handle.outputDir.endswith("a"):
                raise RuntimeError("kill failed")

        monkeypatch.setattr(jobtrack, "cancelGroup", flaky)
        jobtrack.cancelAllGroups(reg)
        assert len(calls) == 2
        assert reg.active_groups() == []
        assert not (tmp_path / "a" / "running.pid").exists()
