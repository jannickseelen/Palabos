#!/bin/bash

set -e

cd ~/Palabos

git fetch --all

git reset --hard origin/master
