#!/usr/bin/env bash

set -xe

mkdir -p include
cd include
git clone https://github.com/lvandeve/lodepng.git && mv lodepng/lodepng.cpp lodepng/lodepng.c
cd ..

CFLAGS="-std=c99 -c -O3 -march=native -Wall -Wextra -Wpedantic -Wshadow -Wfloat-conversion -Wno-strict-aliasing" 

$CC $CFLAGS -o lodepng.o include/lodepng/lodepng.c
$CC $CFLAGS filters.c
$CC $CFLAGS resize.c
$CC $CFLAGS helper.c
$CC $CFLAGS example.c

$CC -o resize *.o

mkdir -p build
mv resize build

rm *.o
