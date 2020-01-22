#!/usr/bin/env bash

podman run -w /root --rm -e DISABLE_AUTH=true \
       --mount type=bind,src="$PWD",dst=/root \
       docker.io/eamon/2019transients:v20200102 \
       Rscript -e "knitr::knit(\"reactivity-simulation-tidy.Rmd\")"
