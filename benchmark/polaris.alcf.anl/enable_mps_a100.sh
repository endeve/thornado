#!/bin/bash -l

export CUDA_MPS_PIPE_DIRECTORY=/tmp/nvidia-mps
export CUDA_MPS_LOG_DIRECTORY=/tmp/nvidia-log

CUDA_VISIBLE_DEVICES=0 nvidia-cuda-mps-control -d
echo "start_server -uid $( id -u )" | nvidia-cuda-mps-control
