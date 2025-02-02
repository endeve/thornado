#!/usr/bin/env bash

sudo apt update
sudo apt -y upgrade
sudo apt -y install \
  libhdf5-dev \
  libopenmpi-dev \
  liblapack-dev \
  libblas-dev
