# transient-indicators

[![DOI](https://zenodo.org/badge/235404395.svg)](https://zenodo.org/badge/latestdoi/235404395)

These files are available to reproduce the results for the manuscript
"Transient indicators of tipping points in infectious diseases." Some
of the results can be reproduced using the Mathematica Notebook, and
the rest can be reproduced using the R markdown notebook. The
Dockerfile can be used to build a Linux container image with a
suitable version of R and all required packages for running the R
markdown notebook. The bash script `build-image.sh` could be used to
generate such an image with podman version 1.6.2 on Ubuntu 19.10. An
already built image has also been pushed to Docker Hub. The bash
script `run-simulations.sh` could be used to run the R markdown
notebook with that image using podman. The directory
`example-output-2020-01-22` contains files produced by running that
script. Producing them took about 75 minutes using eight 2.1GHz cores.
