#!/bin/bash
orterun -bynode -np 1 --hostfile openmpi.config.file /http/R-patched5/bin/R --slave < run-problem-openmpi.R