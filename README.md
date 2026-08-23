# Compositron - Modular Multi-Implementation Analysis Software Compositor for Multi-Technique Positron Annihilation Spectroscopy

## Description
`Compositron` is a compositor of analysis software for positron annihilation spectroscopy.
It composes modules for different techniques, each one implemented multiple times in different programming languages.
The focus of this project is primarily on the technical part, not the physical one.
The modules can be used directly as libraries, or in a GUI for graphical data evaluation.
The core computational libraries and the GUIs are designed to be modular, allowing them to be combined arbitrarily.

## Roadmap

**Short Term**

| Feature              | Rust  | C++   | Julia | Python | C     |
| :---                 | :---: | :---: | :---: | :---:  | :---: |
| SLOPE importer       | [x]   | [x]   | [x]   | [x]    | [x]   |
| DBS analysis         | [x]   | [x]   | [x]   | [x]    | [ ]   |
| MePS importer        | [x]   | [x]   | [x]   | [ ]    | [ ]   |
| PALS analysis        | [x]   | [x]   | [x]   | [ ]    | [ ]   |
| CDBS analysis        | [x]   | [x]   | [ ]   | [x]    | [ ]   |
| Measurement Handling | [ ]   | [ ]   | [ ]   | [ ]    | [ ]   |

**Long Term (in only one of the above languages)**

- Profile analysis (VEPFIT, LIMPID)
- Advanced PALS analysis (MELT, CONTIN)
- Advanced CDBS analysis (In-Flight)
- GPU acceleration
- Listmode processing tools
- Waveform processing tools
- AMOC analysis
- ACAR analysis

## Architectures

Currently, the project is being developed in the following programming languages:

- Python
- C++
- Julia
- Rust
- C

Since core library and GUI are not generally of the same language, an extra interface layer is needed for interoperability.
There are different ways to do that, and (if time allows) they will all be implemented in this project to benchmark them against each other:

- Making a server process that can use the core library directly and communicate with the GUI via a common IPC protocol, such as HTTP or [gRPC](https://grpc.io/).
- Putting foreign language bindings on top of the libraries. For example, Python bindings on the C library allow a Python GUI to use the C library directly.

Another approach is to write monoliths, i.e., GUI and engine are of the same language, run in (possibly) just a single process and are closely integrated.

Currently, there is only one GUI: JavaScript (using [Angular](https://angular.dev/) and [Electron](https://www.electronjs.org/)) and it communicates with the cores via HTTP.

The implementation of GUIs with the following frameworks is planned as well:

- Qt (C++)
- PyQt (Python)

GPU acceleration is planned at least for the inversion of lifetime spectra and at least in CUDA.

To use any of the components, refer to the READMEs in the respective folders.

## Benchmarks

Coming soon...

## Disclaimers

Each part of this project is in the development stage and not guaranteed to work or produce accurate results. Use at your own peril.

Changelogs are nice, but currently not used in this project. They will be, if anything is mature enough to be published or if enough colleagues start using the GUIs offered by this project.
