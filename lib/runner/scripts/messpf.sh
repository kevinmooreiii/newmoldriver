#!/usr/bin/env bash

export OMP_NUM_THREADS=10
messpf pf.inp pf.out >> stdout.log &> stderr.log
