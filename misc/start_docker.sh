#!/usr/bin/env bash


nvidia-docker run --rm --ipc=host -p 8888:8888 \
  -v ~/Projects/Neonympha_classification:/workspace \
  butterflyology/fastai_jupyter jupyter notebook
