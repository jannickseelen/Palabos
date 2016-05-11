#!/bin/bash

set -e

cd ~/Palabos/palabos-v1.5r1/

gitpull.sh
compile.sh
run.sh
