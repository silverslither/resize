#!/usr/bin/env bash

set -x

mkdir -p include
cd include
test -f wuffs-v0.4.c || curl -O https://raw.githubusercontent.com/google/wuffs/refs/heads/main/release/c/wuffs-v0.4.c
test -f fpng.cpp || curl -O https://raw.githubusercontent.com/richgel999/fpng/refs/heads/main/src/fpng.cpp
test -f fpng.h || curl -O https://raw.githubusercontent.com/richgel999/fpng/refs/heads/main/src/fpng.h
cd ..

set -e

CFLAGS="-std=c99 -c -O3 -march=native -fno-strict-aliasing -Wall -Wextra -Wpedantic -Wshadow -Wfloat-conversion"
CXXFLAGS="-std=c++11 -c -O3 -march=native -fno-strict-aliasing -Wall -Wextra -Wpedantic -Wshadow -Wfloat-conversion -Wno-type-limits -Wno-tautological-constant-out-of-range-compare"

$CC $CFLAGS filters.c
$CC $CFLAGS resize.c
$CC $CFLAGS colour.c
$CXX $CXXFLAGS example.cpp

$CXX -lm -o resize *.o -s

mkdir -p build
mv resize build

rm *.o
