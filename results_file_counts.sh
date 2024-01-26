#!/usr/bin/env bash
for x in results/default results/original results/short results/blastn*; do echo $x $(find $x -name "*.csv"| wc -l); done
