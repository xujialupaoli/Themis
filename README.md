# Themis
<<<<<<< HEAD
Themis is a fast and robust metagenomic profiler that achieves high accuracy across ultra-low to high sequencing depths.
=======

**Themis** is a metagenomic profiling tool that integrates a customized version of **ganon** and **ggcat** with a reproducible, automatic post-processing pipeline.

It is designed to look and behave like a standalone, installable software package:

- One **conda** package includes:
  - patched `ganon` (for reference database build and initial classification)
  - `ggcat` (for compact colored DBG indexing and querying)
  - Themis Python pipeline (for filtering, sub-database construction, length correction, adaptive weighting, and final taxonomic profiles)
- One command to build databases.
- One command to generate species-level abundance profiles from raw reads.

Themis is released under the **GPL-3.0** license, and is compatible with ganon and ggcat licenses.

---

## Features

- **Integrated pipeline**
  - `ganon build-custom` to build reference databases
  - `ganon classify + report` to obtain taxonomic profiles
  - Adaptive reference subset selection based on ganon predictions
  - `ggcat build + query` on the filtered reference
  - Length-corrected abundance estimation
  - Hybrid prediction: DBG-based vs ganon-based, with automatic weight selection

- **Production-style CLI**
  - `themis-build-custom`: wrapped ganon database builder
  - `themis-profile`: run the full profiling workflow for paired-end or single-end data

- **Reproducible and self-contained**
  - All components (ganon, ggcat, Python scripts) are bundled in a single conda package.
  - No manual symlink hacks; all intermediate files are managed systematically.

---

## Installation

> Note: Replace `<your-conda-channel>` with your actual channel name once you upload the conda package.

```bash
mamba install -c <your-conda-channel> themis
>>>>>>> 057bdf5 (Themis v0.1.0: initial public release)
