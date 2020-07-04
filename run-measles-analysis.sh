#!/usr/bin/env bash

podman run -w /root --rm -e DISABLE_AUTH=true \
       --mount type=bind,src="$PWD",dst=/root \
       docker.io/eamon/2019transients:v20200626 \
       Rscript -e "knitr::knit(\"measles-outbreaks.Rmd\")"
