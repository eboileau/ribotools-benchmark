# Ribotools benchmark

A workflow to benchmark tools for differential translation efficiency (TE) analysis and a companion repository to

> Rp-Bp and the Ribotools box for translatome profiling, Etienne Boileau, Maja Bencun, Francesca Tuorto, Isabel S. Naarmann-de Vries, Pauline Fahjen, Philipp Mertins, and Christoph Dieterich.

## Scope and documentation

This repository is a supplementary companion, not a standalone package. It's key use case is to benchmark selected TE software tools using the simulated data from [Oertlin et al.](https://academic.oup.com/nar/article/47/12/e70/5423604) and a subset of the data from [Boileau et al.](https://zenodo.org/records/17899627). Although this is a well-documented space with plenty of published benchmarks, as of 2026, many tools have been deprecated and are now obsolete or unmaintained. A new wave of E2E software such as [Ribotools](https://ribotools.readthedocs.io) are now making the analysis of TE accessible, reproducible, and interpretable. For more details, results and discussion, consult the manuscript.

## Installation

```bash
# install snakemake into an environment named "benchmark"
mamba create -c conda-forge -c bioconda -c nodefaults -n benchmark snakemake-minimal
# install conda into this environment
mamba activate benchmark
mamba install -c conda-forge conda
```

**Note:** Setting `conda config --set channel_priority strict` may break environment specs for selected tools. Use the default settings or
explicitly set channel priority to `flexible`.

## Usage

### Configuration

To configure this workflow, modify the following files:

+ [config/config.yaml](config/config.yaml): workflow-specific parameters

This use case is self-documenting. If you want to add a new tool to the benchmark, update the configuration file, add an environment file *e.g.* [workflow/envs/ribotools.yaml](workflow/envs/ribotools.yaml) and a corresponding script *e.g.* [workflow/scripts/run_ribotools](workflow/scripts/run_ribotools.R). More complex use cases may require additional modifications. Community contributions are welcome.

#### Profile

Snakemake scheduling hints are set in a default profile under [workflow/profiles/default/config.yaml](workflow/profiles/default/config.yaml). In addition, processes memory and CPU usage (*e.g.* BLAS and OMP used via parallel/BiocParallel by some R packages) are capped at the OS level in Snakemake rules.

### How to run the workflow

```bash
snakemake --profile workflow/profiles/default
```

### Output
