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
