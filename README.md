# Themis

**Themis** is a robust and accurate species-level metagenomic profiler.

[![BioConda Install](https://img.shields.io/conda/dn/bioconda/themis.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/themis)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/themis/badges/version.svg)](https://anaconda.org/bioconda/themis)
[![License](https://img.shields.io/github/license/xujialupaoli/Themis)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/xujialupaoli/Themis)](https://github.com/xujialupaoli/Themis/releases)
[![GitHub Downloads](https://img.shields.io/github/downloads/xujialupaoli/Themis/total.svg?style=social&logo=github&label=Download)](https://github.com/xujialupaoli/Themis/releases)


## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick start](#quick-start)
  - [1. Build custom reference database](#1-build-custom-reference-database)
  - [2. Run profiling](#2-run-profiling)
- [Input files](#input-files)
- [Output files](#output-files)
- [Command reference](#command-reference)
- [Tips and known issues](#tips-and-known-issues)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## Features

- **Commands**
  - `themis build-custom` Build custom themis databases.
  - `themis profile` Profile reads against custom databases.




---

## Installation

```bash
mamba install -c bioconda themis
