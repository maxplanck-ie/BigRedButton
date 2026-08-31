# Changelog

## [0.8.1](https://github.com/maxplanck-ie/BigRedButton/compare/v0.8.0...v0.8.1) (2026-08-31)


### Bug Fixes

* route Aviti samba QC copies to their machine's own facility folder ([#149](https://github.com/maxplanck-ie/BigRedButton/issues/149)) ([382f5bc](https://github.com/maxplanck-ie/BigRedButton/commit/382f5bce07aaf229a8d5fb08730acf6fa7767420))

## [0.8.0](https://github.com/maxplanck-ie/BigRedButton/compare/v0.7.1...v0.8.0) (2026-08-31)


### Features

* Dispatch a flowcell's library-groups concurrently through a bounded thread pool ([#145](https://github.com/maxplanck-ie/BigRedButton/issues/145)) ([061c90a](https://github.com/maxplanck-ie/BigRedButton/commit/061c90a7bb8169caba8161f5653dfd24bdbf431b))

## [0.7.1](https://github.com/maxplanck-ie/BigRedButton/compare/v0.7.0...v0.7.1) (2026-08-24)


### Bug Fixes

* nest aviti scanning, output, and logs under the matching serial-ID subdir ([#143](https://github.com/maxplanck-ie/BigRedButton/issues/143)) ([827bb51](https://github.com/maxplanck-ie/BigRedButton/commit/827bb51db3759f0bdefda7804f20d952968f1f65))

## [0.7.0](https://github.com/maxplanck-ie/BigRedButton/compare/v0.6.0...v0.7.0) (2026-08-12)


### Features

* process external RELACS ChIP-Seq and fix RELACS QC-sharing bugs ([#137](https://github.com/maxplanck-ie/BigRedButton/issues/137)) ([9a9d932](https://github.com/maxplanck-ie/BigRedButton/commit/9a9d932ac2de721062b43cd6b9bf39b3ced8ded3))

## [0.6.0](https://github.com/maxplanck-ie/BigRedButton/compare/v0.5.0...v0.6.0) (2026-08-11)


### Features

* add install_brb.sh for versioned conda env installs from tags ([#138](https://github.com/maxplanck-ie/BigRedButton/issues/138)) ([05ca608](https://github.com/maxplanck-ie/BigRedButton/commit/05ca6082a0b44fb7cbf9a2e4cf256bef5aa47c9f))

## [0.5.0](https://github.com/maxplanck-ie/BigRedButton/compare/v0.4.1...v0.5.0) (2026-08-10)


### Features

* Add -s/--sequencer to restrict polling to one platform's baseData/logPath ([#134](https://github.com/maxplanck-ie/BigRedButton/issues/134)) ([0782d96](https://github.com/maxplanck-ie/BigRedButton/commit/0782d96ade6faaecbea3b15d8581419a5f5ba1b0))

## [0.4.1](https://github.com/maxplanck-ie/BigRedButton/compare/v0.4.0...v0.4.1) (2026-07-23)


### Bug Fixes

* Scope Aviti detection glob to the current run directory ([#129](https://github.com/maxplanck-ie/BigRedButton/issues/129)) ([e164d4d](https://github.com/maxplanck-ie/BigRedButton/commit/e164d4d6a1bf08329414d99500d7acc03371a721))

## [0.4.0](https://github.com/maxplanck-ie/BigRedButton/compare/v0.3.0...v0.4.0) (2026-07-06)


### Features

* built-in Levenshtein distance ([#123](https://github.com/maxplanck-ie/BigRedButton/issues/123)) ([46b4ea3](https://github.com/maxplanck-ie/BigRedButton/commit/46b4ea3572cfa84938551fe74b053b6182573929))

## [0.3.0](https://github.com/maxplanck-ie/BigRedButton/compare/v0.2.0...v0.3.0) (2026-06-10)


### Features

* duplication rate 0 for samples irretrievable due to poor naming ([82e7498](https://github.com/maxplanck-ie/BigRedButton/commit/82e749898f7e5d643434f444e44c57e3faa8ccc1))
* duplication rate 0 for samples irretrievable due to poor naming ([#120](https://github.com/maxplanck-ie/BigRedButton/issues/120)) ([82e7498](https://github.com/maxplanck-ie/BigRedButton/commit/82e749898f7e5d643434f444e44c57e3faa8ccc1))

## [0.2.0](https://github.com/maxplanck-ie/BigRedButton/compare/v0.1.0...v0.2.0) (2025-10-27)


### Features

* Added a function to detect the sequencing type (decision maker for subsequent sequencing type check. ([132ce0c](https://github.com/maxplanck-ie/BigRedButton/commit/132ce0c73905efd1ae8423812177c02c0c6045db))
* aviti integration and SAMBA structure update ([#114](https://github.com/maxplanck-ie/BigRedButton/issues/114)) ([132ce0c](https://github.com/maxplanck-ie/BigRedButton/commit/132ce0c73905efd1ae8423812177c02c0c6045db))
* Implemented Aviti flowcell ID recognition and Parkour communication for downstream analysis. ([132ce0c](https://github.com/maxplanck-ie/BigRedButton/commit/132ce0c73905efd1ae8423812177c02c0c6045db))
* Introduced new SAMBA structure for improving data management. ([132ce0c](https://github.com/maxplanck-ie/BigRedButton/commit/132ce0c73905efd1ae8423812177c02c0c6045db))


### Bug Fixes

* Added --snakemakePath option to enable the 10X_snakePipe workflow. ([132ce0c](https://github.com/maxplanck-ie/BigRedButton/commit/132ce0c73905efd1ae8423812177c02c0c6045db))
* Updated adapters for the RELACS pipeline. ([132ce0c](https://github.com/maxplanck-ie/BigRedButton/commit/132ce0c73905efd1ae8423812177c02c0c6045db))

## [0.1.0](https://github.com/maxplanck-ie/BigRedButton/compare/v0.0.15...v0.1.0) (2025-10-06)


### Features

* release please ([#112](https://github.com/maxplanck-ie/BigRedButton/issues/112)) ([21887fe](https://github.com/maxplanck-ie/BigRedButton/commit/21887fed96807f63c41748ded2ffc16ff4b14dd4))


### Bug Fixes

* drop reuse of quote character in f-strings for overall python compatibility ([21887fe](https://github.com/maxplanck-ie/BigRedButton/commit/21887fed96807f63c41748ded2ffc16ff4b14dd4))
