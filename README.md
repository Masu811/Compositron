# Compositron

## Description

`Compositron` is a High-Performance All-In-One Data Analysis Toolkit for Multi-Technique Positron Annihilation Spectroscopy.
It is written in Rust, but currently also contains some modules additionally in Python, Julia and C++ for performance comparison.
The proclaimed goal of this project is high performance at all costs, which makes the source code harder to maintain. The current plan is that `Compositron` will be retired, if PANDORA (to be published) and NEPOMUC's Python package suite ever reach the same level of performance and sophistication.
Until then, `Compositron` continues growing.
A Tauri GUI is definitely planned.

## Roadmap

| Feature                      | Rust  | C++   | Julia | Python |
| :---                         | :---: | :---: | :---: | :---:  |
| DBS analysis                 | ✓     | ✓     | ✓     | ✓      |
| Basic CDBS analysis          | ✓     | (✓)   | (✓)   | (✓)    |
| Advanced CDBS analysis       |       |       |       |        |
| Basic PALS analysis          | ✓     | ✓     | ✓     | ✓      |
| Advanced PALS analysis       |       |       |       |        |
| ACAR analysis                |       |       |       |        |
| AMOC analysis                |       |       |       |        |
| Listmode analysis            |       |       |       |        |
| Waveform analysis            | ✓     |       |       |        |
| Profile analysis (Diffusion) |       |       |       |        |
| Synthetic data generator     |       |       |       |        |
| Multi-Threading              |       |       |       |        |
| GPU acceleration             |       |       |       |        |
| SIMD vectorization           |       |       |       |        |

## Benchmarks

Coming soon...

## Disclaimers

Each part of this project is in the development stage and not guaranteed to work or produce accurate results. Use at your own peril.

Changelogs are nice, but currently not used in this project. They will be, if anything is mature enough to be published or if enough colleagues start using the GUIs offered by this project.
