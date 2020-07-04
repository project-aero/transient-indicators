FROM docker.io/rocker/tidyverse:3.6.1
MAINTAINER Eamon O'Dea <[last name without apostrophe]35@gmail.com>

RUN install2.r --skipinstalled --error cowplot \
ggforce \
padr \
pomp \
spaero \
tictoc \
zoo
