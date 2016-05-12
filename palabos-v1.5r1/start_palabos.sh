#!/bin/bash

set -e

chmod +x clean.sh
./clean.sh
chmod +x gitpull.sh
./gitpull.sh
chmod +x compile.sh
./compile.sh
chmod +x run.sh
./run.sh
