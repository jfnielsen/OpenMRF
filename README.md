![OpenMRF](OpenMRF_banner.png)

# Documentation

**[OpenMRF: Official Documentation & Blog](https://tomgriesler.github.io/openmrf_documentation/)**

# Introduction

OpenMRF is a modular and vendor-neutral framework for Magnetic Resonance Fingerprinting (MRF) built on the open-source [Pulseq](https://pulseq.github.io) standard. It is built upon the MATLAB version of Pulseq by Kelvin J. Layton and Maxim Zaitsev ([doi:10.1002/mrm.26235](https://doi.org/10.1002/mrm.26235)). OpenMRF unifies all core components of the MRF workflow within a single MATLAB-based toolbox: flexible sequence generation, automated Bloch-based dictionary simulation, and low-rank image reconstruction. The provided tools support a wide range of contrast preparations and readouts (e.g., spiral, radial, rosette) and include integrated solutions for trajectory calibration, spin-lock modeling, slice profile simulation, and metadata storage. Designed for reproducibility and portability, OpenMRF enables easy deployment of MRF protocols across multiple scanner platforms, including Siemens, GE and United Imaging systems.

# Contents

- `include_miitt/`: Contains the low-rank reconstruction code provided by the MIITT group and Jeffrey Fessler's [MIRT toolbox](https://web.eecs.umich.edu/~fessler/code/). Includes an installation script. **Do not** add this folder manually to your MATLAB path; use the `install_OpenMRF.m` script.
- `include_misc/`: Miscellaneous utilities and helper functions.
- `include_pre_sim_library/`: Library of SLR optimized RF pulse waveforms (including `sigpy`-generated pulses) and flip angle patterns for MRF. Also used to store pre-simulated slice profiles, adiabatic efficiencies and compressed dictionaries.
- `include_pulseq_original/`: Copy of the official Pulseq repo ([GitHub link](https://github.com/pulseq/pulseq), v1.5, 01.04.2025). Includes minor modifications to the plotting functions for visualizing trigger inputs/outputs.
- `include_pulseq_toolbox/`: Contains standard imaging readouts (cartesian, radial, spiral, rosette) combined with various preparation modules (inversion recovery, saturation, MLEV-T2, spin-lock, adiabatic spin-lock, CEST). Also includes simulation tools for MRF dictionary generation.
- `main_sequences/`: Example Pulseq sequences and reconstruction scripts.
- `projects/`: Collection of projects, which were published or presented on conferences or which are currently work in progress.
- `user_specifications/`: User specific definitions (automatically generated via install_OpenMRF.m), MRI system specifications (create a .csv file for your system's gradient limits and timings) and optionally your python cmd specifications for SIGPY pulses.

# System Requirements

- **MATLAB** tested with R2024b on Win11 and Ubuntu 22.04
- **Python** with the `sigpy` package — required for designing SLR and adiabatic RF pulses (e.g., BIR-4)

# Note

This repository includes third-party software distributed under their respective licenses. Please consult the [NOTICE](./NOTICE) file before use. Not all code is MIT-licensed. Users intending commercial use must review the [NOTICE](./NOTICE) file and seek appropriate permissions from the original authors.


_Maximilian Gram: 13.11.2025_
