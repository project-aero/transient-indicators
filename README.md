# transient-indicators

[![DOI](https://zenodo.org/badge/235404395.svg)](https://zenodo.org/badge/latestdoi/235404395)

These files are available to reproduce the results for the manuscript
"Transient indicators of tipping points in infectious diseases." Some
of the results can be reproduced using the Mathematica Notebook, and
the rest can be reproduced using the R markdown notebooks. The
Dockerfile can be used to build a Linux container image with a
suitable version of R and all required packages for running the R
markdown notebook. The bash script `build-image.sh` could be used to
generate such an image with podman version 1.9.3 on Ubuntu 18.04. An
already built image has also been pushed to Docker Hub. The bash
scripts `run-simulations.sh`, `run-measles-analysis.sh`, and `run-pertussis-analysis.sh` could be used to run the R markdown
notebook with that image using podman. The directory
`example-output-2020-06-26` contains files produced by running those
script. Producing them took about 90 minutes using eight 3.47GHz cores.
