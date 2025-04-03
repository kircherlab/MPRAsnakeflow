# Changelog

## [0.4.7](https://github.com/kircherlab/MPRAsnakeflow/compare/v0.4.6...v0.4.7) (2025-04-03)


### Bug Fixes

* assignment barcode counting ([9584802](https://github.com/kircherlab/MPRAsnakeflow/commit/95848023860f9b43e9c5a2e594c0fb9d3eda4205))
* merge.label.py does not return pseudocounts in barcode file ([#187](https://github.com/kircherlab/MPRAsnakeflow/issues/187)) ([023110a](https://github.com/kircherlab/MPRAsnakeflow/commit/023110a5368704e0efb464dbaef1a1d5fdcf6c86))
* report also barcodes that are lower than the min barcodes. ([90f1b14](https://github.com/kircherlab/MPRAsnakeflow/commit/90f1b14d04737723fa5852c809c44e497b2e2d3b)), closes [#185](https://github.com/kircherlab/MPRAsnakeflow/issues/185)

## [0.4.6](https://github.com/kircherlab/MPRAsnakeflow/compare/v0.4.5...v0.4.6) (2025-02-26)


### Bug Fixes

* revert v0.4.5. The fix was already implemented within the common.smk file ([7997440](https://github.com/kircherlab/MPRAsnakeflow/commit/799744050b540d96545ac72bc2280e8e37a2b61e))

## [0.4.5](https://github.com/kircherlab/MPRAsnakeflow/compare/v0.4.4...v0.4.5) (2025-02-24)


### Bug Fixes

* change the length of the sequence we want to check in the attach_sequence case ([#177](https://github.com/kircherlab/MPRAsnakeflow/issues/177)) ([e3db5cf](https://github.com/kircherlab/MPRAsnakeflow/commit/e3db5cf487423126ada0fda5305f42a0b38f25cb))

## [0.4.4](https://github.com/kircherlab/MPRAsnakeflow/compare/v0.4.3...v0.4.4) (2025-02-05)


### Bug Fixes

* qc report count typo fix ([#174](https://github.com/kircherlab/MPRAsnakeflow/issues/174)) ([d78d339](https://github.com/kircherlab/MPRAsnakeflow/commit/d78d3391700c092c90742d1e47dc7b65b4d6d808))

## [0.4.3](https://github.com/kircherlab/MPRAsnakeflow/compare/v0.4.2...v0.4.3) (2025-01-29)


### Bug Fixes

* swap activity threshold and all output file paths in qc_report rule ([#172](https://github.com/kircherlab/MPRAsnakeflow/issues/172)) ([ff1b11c](https://github.com/kircherlab/MPRAsnakeflow/commit/ff1b11cc3baf84a0a0c45201d93b36d8f630d42b))

## [0.4.2](https://github.com/kircherlab/MPRAsnakeflow/compare/v0.4.1...v0.4.2) (2025-01-07)


### Bug Fixes

* swapped barcode output files ([a851ad2](https://github.com/kircherlab/MPRAsnakeflow/commit/a851ad2d74aad0834b853eb7432da10f369f8e21))

## [0.4.1](https://github.com/kircherlab/MPRAsnakeflow/compare/v0.4.0...v0.4.1) (2025-01-06)


### Bug Fixes

* update configs to 0.4 ([bd79b6b](https://github.com/kircherlab/MPRAsnakeflow/commit/bd79b6b22118749a5eec89ec95815b243318dd25))

## [0.4.0](https://github.com/kircherlab/MPRAsnakeflow/compare/v0.3.1...v0.4.0) (2025-01-02)


### ⚠ BREAKING CHANGES

* Major version change because of changes in config that are not backward compatible.

### Features

* Add performance tweaks for resource optimization in workflow rules ([5ed1ef9](https://github.com/kircherlab/MPRAsnakeflow/commit/5ed1ef9fb92acc56c999f26183e5ef0d4caf5be0))
* Allowing also only FW reads with a UMI ([56b2254](https://github.com/kircherlab/MPRAsnakeflow/commit/56b225405bba024779edba11a332872b01ac02de))
* Strand sensitive option ([f315aab](https://github.com/kircherlab/MPRAsnakeflow/commit/f315aabd24ebe6f524b9fb3b008f5af96987775d))
* One DNA or RNA count file across multiple replicates ([ccbacee](https://github.com/kircherlab/MPRAsnakeflow/commit/ccbacee28bd551f9508e80d4e966188f970ba077))

## [0.3.1](https://github.com/kircherlab/MPRAsnakeflow/compare/v0.3.0...v0.3.1) (2024-12-17)


### Bug Fixes

* Wrong experiment count plots in QC report ([#149](https://github.com/kircherlab/MPRAsnakeflow/issues/149)) ([d2be468](https://github.com/kircherlab/MPRAsnakeflow/commit/d2be46891650ff9aaab61f750a4b3bc3b65e3e88))

## [0.3.0](https://github.com/kircherlab/MPRAsnakeflow/compare/v0.2.0...v0.3.0) (2024-11-20)


### ⚠ BREAKING CHANGES

* versioned config ([#140](https://github.com/kircherlab/MPRAsnakeflow/issues/140))

### Features

* versioned config ([#140](https://github.com/kircherlab/MPRAsnakeflow/issues/140))
* MAD outlier removal is completely removed ([#140](https://github.com/kircherlab/MPRAsnakeflow/issues/140))
* default is NO outlier detection (none is not present anymore) ([#140](https://github.com/kircherlab/MPRAsnakeflow/issues/140))
* global config is removed. splits moved now withing mapping in assignment ([#140](https://github.com/kircherlab/MPRAsnakeflow/issues/140))

## [0.2.0](https://github.com/kircherlab/MPRAsnakeflow/compare/v0.1.1...v0.2.0) (2024-11-05)

### ⚠ BREAKING CHANGES

* Support only snakemake >=8.24.1 ([#130](https://github.com/kircherlab/MPRAsnakeflow/pull/130))
* File output formats and locations changed
* Normalization changed which may result in different outputs

### Features
 
 * outlier removal methods ([#132](https://github.com/kircherlab/MPRAsnakeflow/pull/132))
 * No min max length for bbmap. default mapq is 30. ([#131](https://github.com/kircherlab/MPRAsnakeflow/pull/131))
 * IGVF outputs ([#129](https://github.com/kircherlab/MPRAsnakeflow/pull/129))
 * Documentation improvements


### Bug Fixes

## [0.1.1](https://github.com/kircherlab/MPRAsnakeflow/compare/v0.1.0...v0.1.1) (2024-09-30)

### Bug Fixes

* Detach from anaconda ([#122](https://github.com/kircherlab/MPRAsnakeflow/issues/122)) ([16bcea2](https://github.com/kircherlab/MPRAsnakeflow/commit/16bcea2f04190a5965ad1865cf30f6dd44f1b6a0))
* memory resources for bbmap ([#123](https://github.com/kircherlab/MPRAsnakeflow/issues/123)) ([af93f58](https://github.com/kircherlab/MPRAsnakeflow/commit/af93f588e9387ddf91197f5587d36c3481499b38))

## [0.1.0](https://github.com/kircherlab/MPRAsnakeflow/compare/v0.0.1...v0.1.0) (2024-09-18)

First release of MPRAsnakeflow! 

### Feature highlights

* Multiple assignment mapping strategies (BBMap, exact, bwa)
* Optimized assignment for variants with BBMap
* QC report for assignment and experiment workflow
* Barcode count output
* Snakemake 8 support
* Extended documentation: https://mprasnakeflow.readthedocs.io


## older development


### ⚠ BREAKING CHANGES

* latest development for new release ([#133](https://github.com/kircherlab/MPRAsnakeflow/issues/133))
* pseudocounts where not used correctly when RNA or DNA set to 0
* DNA and RNA join correction

### Features

* Add assignment_merge thread configuration ([26e68c2](https://github.com/kircherlab/MPRAsnakeflow/commit/26e68c26f315c524cf28692d636127fbf3bdeb2b))
* better assignment BC statistics ([00187e6](https://github.com/kircherlab/MPRAsnakeflow/commit/00187e689b2fad10fd317aa2efbd0214fad14434))
* configurable min mapping quality ([28045ae](https://github.com/kircherlab/MPRAsnakeflow/commit/28045aea23d6fa03f3883b3dc44b3cbc3e8f6205))
* extending figure width ([8bf81c4](https://github.com/kircherlab/MPRAsnakeflow/commit/8bf81c45e45f9b4c23856c0915bd527f9699b6cd))
* faster design check ([315b402](https://github.com/kircherlab/MPRAsnakeflow/commit/315b402499d92850382d4110e153602020381e8a))
* fastq-join implementation ([aaf5315](https://github.com/kircherlab/MPRAsnakeflow/commit/aaf5315364ebb3e3117c3996c2fc357aa9c4d595))
* latest development for new release ([#133](https://github.com/kircherlab/MPRAsnakeflow/issues/133)) ([bdfc557](https://github.com/kircherlab/MPRAsnakeflow/commit/bdfc557a64cecc19d1d86eead8bdb691a1ff2166))
* make filtering consistent ([5f7a4c5](https://github.com/kircherlab/MPRAsnakeflow/commit/5f7a4c5a2a3389a75b8d6b7e9aaf34485127b3a4))
* master variant table ([6bda47c](https://github.com/kircherlab/MPRAsnakeflow/commit/6bda47c78021bc1728bb81a716f5e6daaf6ac084))
* new final output file with merged replicates ([66cf017](https://github.com/kircherlab/MPRAsnakeflow/commit/66cf0172cb6b556e507be4daabf7e859447787f3))
* only link assignment fasta when possible ([d7d3822](https://github.com/kircherlab/MPRAsnakeflow/commit/d7d3822933c98d790f3c96bcbfdef1a7ea70c7df)), closes [#50](https://github.com/kircherlab/MPRAsnakeflow/issues/50)
* remove space, speedup BC extraction ([70e9bd0](https://github.com/kircherlab/MPRAsnakeflow/commit/70e9bd06b91ccb37333e0a69c47917a5eacbf639))
* replace merging by NGmerge ([0aa8cad](https://github.com/kircherlab/MPRAsnakeflow/commit/0aa8cad6884a953f9c89a2fdd7af397e4e9ccf3e))
* snakemake 8 compatibility ([cf38ed9](https://github.com/kircherlab/MPRAsnakeflow/commit/cf38ed9de68367d0d1700ccff262e91ad6f1fbc0))
* snakemake 8 ready with workflow profile ([d637e1f](https://github.com/kircherlab/MPRAsnakeflow/commit/d637e1fdbebfca0616d944101898fbf522df9c82))
* statistic for assignment workflow ([10c3b26](https://github.com/kircherlab/MPRAsnakeflow/commit/10c3b2677ada59925ddd3de777f7488c9a20e981))
* using reverese compelment BCs ([d009a6c](https://github.com/kircherlab/MPRAsnakeflow/commit/d009a6c3de7de50a210479b73f5d41969287e234))


### Bug Fixes

* batch size issue in sort ([487ba8c](https://github.com/kircherlab/MPRAsnakeflow/commit/487ba8ce059517030fcab3708c3cea40ac210f7e))
* correct use of assignment configs ([58b64f1](https://github.com/kircherlab/MPRAsnakeflow/commit/58b64f1e753477f7410233ac546701ddbd60f9f2))
* corrected qc_report_assoc ([afb0127](https://github.com/kircherlab/MPRAsnakeflow/commit/afb012750bc1c3c39f2348b283c23ff97695f672))
* Detach from anaconda ([#122](https://github.com/kircherlab/MPRAsnakeflow/issues/122)) ([16bcea2](https://github.com/kircherlab/MPRAsnakeflow/commit/16bcea2f04190a5965ad1865cf30f6dd44f1b6a0))
* DNA and RNA join correction ([7214743](https://github.com/kircherlab/MPRAsnakeflow/commit/7214743008dc6796077e45e62646174ffaf52290))
* filter config ([38ee37e](https://github.com/kircherlab/MPRAsnakeflow/commit/38ee37ecfcf4a71b840575504811512e0d64609a))
* issue with stats and asisgnment ([d935fa1](https://github.com/kircherlab/MPRAsnakeflow/commit/d935fa1f62825dfdcd2cd77e4c73bc37686519a0))
* memory resources for bbmap ([#123](https://github.com/kircherlab/MPRAsnakeflow/issues/123)) ([af93f58](https://github.com/kircherlab/MPRAsnakeflow/commit/af93f588e9387ddf91197f5587d36c3481499b38))
* plots per insert only used last experiment. not all. ([c2fd82b](https://github.com/kircherlab/MPRAsnakeflow/commit/c2fd82b6d4b545cc3a1acc5ecb145eb3c93af49d))
* pseudocounts where not used correctly when RNA or DNA set to 0 ([d2483f9](https://github.com/kircherlab/MPRAsnakeflow/commit/d2483f9c7724e0b63cec4f251519d449831ecf04))
* remove illegal characters from reference ([0ebee81](https://github.com/kircherlab/MPRAsnakeflow/commit/0ebee81d74f3f6170ce4b8083e18c746550154db))
* rename barcoe output header ([635f043](https://github.com/kircherlab/MPRAsnakeflow/commit/635f0431c78d3d5bf9b77a16f6ce26d9ff6c82c2))
* rule make_master_tables fix ([df42845](https://github.com/kircherlab/MPRAsnakeflow/commit/df42845b6dfa9a7b64f187b38f1f15518f3e4a31))
* statistic total counts ([6381b92](https://github.com/kircherlab/MPRAsnakeflow/commit/6381b928fd6c14eb16801a459b8546fa37004c74))
* typo in report ([ace8cca](https://github.com/kircherlab/MPRAsnakeflow/commit/ace8ccacb3d7ece04af43c9b0b1dc9c9c087a2c4))
* upgrade code to new pandas version ([aaea236](https://github.com/kircherlab/MPRAsnakeflow/commit/aaea236bc83f459e7a6c2d3fee96d49c79762325))
* using correct threads ([6dcad7d](https://github.com/kircherlab/MPRAsnakeflow/commit/6dcad7d34173f37d4538644b1ba0d918afd8f149))
* using multiple fastq inputs in counts ([95935cf](https://github.com/kircherlab/MPRAsnakeflow/commit/95935cfe69956ca50307a9c6a774c4b96dff860f))
