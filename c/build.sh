#!/usr/bin/env bash

set -x

mkdir -p include
cd include
test -f wuffs-v0.4.c || curl -O https://raw.githubusercontent.com/google/wuffs/refs/heads/main/release/c/wuffs-v0.4.c
test -f fpng.cpp || curl -O https://raw.githubusercontent.com/richgel999/fpng/refs/heads/main/src/fpng.cpp
test -f fpng.h || curl -O https://raw.githubusercontent.com/richgel999/fpng/refs/heads/main/src/fpng.h
cd ..

set -e

test -z $CC && CC=clang
test -z $CXX && CXX=clang++

CFLAGS="-std=c99 -c -O3 -march=native -fno-strict-aliasing -Wall -Wextra -Wpedantic -Wshadow -Wfloat-conversion"
CXXFLAGS="-std=c++17 -c -O3 -march=native -fno-strict-aliasing -Wall -Wextra -Wpedantic -Wshadow -Wfloat-conversion"

$CC $CFLAGS filters.c
$CC $CFLAGS resize.c
$CC $CFLAGS colour.c
$CXX $CXXFLAGS -Wno-type-limits -Wno-tautological-constant-out-of-range-compare -Wno-unused-function include/fpng.cpp
$CXX $CXXFLAGS -Wold-style-cast example.cpp

$CXX -lm -o resize *.o

mkdir -p build
mv resize build

rm *.o
