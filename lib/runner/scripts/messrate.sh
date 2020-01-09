#!/usr/bin/env bash

export OMP_NUM_THREADS=10
mess mess.inp rate.out >> stdout.log &> stderr.log
