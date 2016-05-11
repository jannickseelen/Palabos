#!/bin/bash

set -e

cd ~/Palabos/palabos-v1.5r1/

chmod +x gitpull.sh
gitpull.sh
chmod +x compile.sh
compile.sh
chmod +x run.sh
run.sh
