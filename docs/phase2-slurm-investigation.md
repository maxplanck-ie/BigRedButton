# Phase 2 prerequisite investigation: does snakemake clean up its own Slurm jobs?

**Date:** 2026-08-24
**Host:** `rapidus` (as `pipegrp`)
**Answers:** Task 1 of `/tmp/2026-08-24-brb-phase2-job-cancellation-plan.md`

## Environment

- `~/configs/BigRedButton.ini`: `snakemakeWorkflowBaseDir=/localenv/pipegrp/anaconda/miniconda3/envs/snakepipesV3.2.0`
- snakemake **8.11.0**, snakePipes installed under the same conda env (import name `snakePipes`, not `snakepipes`)
- Installed executor plugins (`pip list`): `snakemake-executor-plugin-cluster-generic 1.0.9` only — no `snakemake-executor-plugin-slurm`, confirming the wiki's statement that the built-in Slurm executor is unsupported here.
- Shared profile: `/data/repository/misc/snakemake_profile/config.yaml`
  - `executor: cluster-generic`
  - `cluster-generic-submit-cmd`: wraps `sbatch --partition=... --cpus-per-task=... --mem-per-cpu=... --job-name=... --output=... --parsable` (bare job ID on stdout, not Slurm's default `Submitted batch job <id>` text)
  - `cluster-generic-cancel-cmd: "ccancel.sh"` — a one-line wrapper (`module load slurm && scancel "$@"`) shipped alongside the profile
- `sacct` is not on `PATH` even after `module load slurm` (`sacct: command not found`) — confirms the wiki's "our slurm setup lacks accounting features."
- **Host note:** neither `module` nor the Slurm client binaries (`squeue`/`sbatch`/`scancel`) are available in a plain non-interactive SSH shell — `~/.bashrc` early-returns for non-interactive shells, and module init lives in `/etc/profile.d/modules.sh`, which only non-interactive *login* shells or an explicit `source` pick up. Any future BRB code or ops tooling that shells out to `sbatch`/`scancel` etc. on `rapidus` must explicitly `source /etc/profile.d/modules.sh && module load slurm` first — it cannot rely on the ambient environment `BRB` itself normally runs in.

## Job-ID observability (Step 3)

| Source | Usable? | Notes |
|---|---|---|
| (a) driver stdout | **Yes — the only reliable one** | Exact log line, verbatim, from a live run: `Submitted job 1 with external jobid '10874970'.` (source: `snakemake_executor_plugin_cluster_generic/__init__.py:176`, format string `"Submitted {} {} with external jobid '{}'."`). Because the profile passes `sbatch --parsable`, the captured ID is a bare integer, not `Submitted batch job <id>`. |
| (b) `.snakemake` state | No | No `.snakemake` directories found anywhere under `/rapidus` at rest — `tidyUpABit` (`BRB/PushButton.py:243`) deletes it after a successful run, and it wouldn't exist yet before a run starts. |
| (c) `clusters_logs` filenames | No | None found under `/rapidus`, same lifecycle reason as (b). |
| (d) `sacct` | No | Not on `PATH`, confirmed above. |

**Conclusion for Task 4:** job IDs must come from parsing the driver's own stdout, as the plan assumed. Update the regex/test fixture to the verbatim line above: `Submitted job \d+ with external jobid '(\d+)'\.`

## The decisive experiment (Step 4)

Ran a trivial `sleep 900` Snakefile under the exact shared profile (`snakemake --profile /data/repository/misc/snakemake_profile -j 1`), twice:

**Run 1 — bare SIGTERM to the driver PID:**
- Job `10874970` submitted, seen as `R` (running) in `squeue`.
- `kill -TERM <driver_pid>` sent directly (no process group).
- 20s later: job `10874970` still `R` (running) in `squeue`; driver process **still alive**.
- `driver.log` tail: `Submitted job 1 with external jobid '10874970'.` / `Will exit after finishing currently running jobs (scheduler).` — then nothing further; the driver just sits there.

**Run 2 — driver in its own session (`setsid`), SIGTERM to the process group:**
- Job `10874971` submitted, seen as `R` (running) in `squeue`.
- `kill -TERM -<driver_pid>` sent to the whole process group (what `runManagedSubprocess`'s `start_new_session=True` + group signal would actually do).
- 20s later: job `10874971` still `R` (running) in `squeue`; driver process **still alive**.
- Identical log tail: `"Will exit after finishing currently running jobs (scheduler)."`

Both variants produced the same outcome. Cross-checked against the snakemake 8.11.0 source on this host to explain *why*:

- `Scheduler.__init__` installs `signal.signal(signal.SIGTERM, self.exit_gracefully)` (`scheduler.py:145`).
- `exit_gracefully` only sets `self._user_kill = "graceful"` (`scheduler.py:471-474`) — it does not raise, does not call `executor.cancel()`.
- In the scheduling loop, when `user_kill == "graceful"` and jobs are still `running`, snakemake logs `"Will exit after finishing currently running jobs (scheduler)."` and `continue`s the loop (`scheduler.py:204-217`) — it does **not** shut down or cancel while anything is still running.
- `executor.cancel()` — the only path that invokes `cancel_jobs()` (and therefore the profile's `ccancel.sh` → `scancel`) — is reached solely from a separate `except (KeyboardInterrupt, SystemExit): ... self._executor.cancel()` handler around the top-level scheduling call (`scheduler.py:318-323`). Installing a custom `SIGTERM` handler means Python's default terminate-on-SIGTERM behavior never fires and nothing re-raises `SystemExit` from `exit_gracefully`, so this path is never reached by a `SIGTERM`.

Cleanup: both probe jobs `scancel`'d, both driver processes killed, `~/scratch/brb-sigterm-probe` removed. `squeue -u pipegrp` empty afterward. No other users' jobs (`10874957`, the unrelated `LONGUS2-...` job seen throughout) were touched.

## Conclusion

**Conclusion A confirmed, both by live experiment and by source inspection:** snakemake does **not** cancel its submitted Slurm jobs when its driver receives `SIGTERM` (bare or process-group). It stops submitting new work and waits indefinitely for already-running jobs to finish on their own. A killed/crashed BRB would leave orphaned Slurm jobs running with no live process left to report their results.

**Decision: proceed with Tasks 2-12 of the Phase 2 plan as written.** Phase 2 needs its own job-ID tracking (parsed from driver stdout, per the format above) plus explicit `scancel` in the crash handler — `SIGTERM` delivery to the driver alone is not sufficient.

**Task 4 update:** use the verbatim line format captured above; `parseJobIds`'s regex should be `Submitted job \d+ with external jobid '(\d+)'\.` (anchor on the quoted group only, since the bare-integer form is specific to `--parsable`).
