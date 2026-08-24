import json
import os
import threading

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
