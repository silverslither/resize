#!/usr/bin/env bash

set -xe

mkdir -p include
cd include
git clone https://github.com/lvandeve/lodepng.git && mv lodepng/lodepng.cpp lodepng/lodepng.c
cd ..

CFLAGS="-std=c99 -c -O3 -mfma -mavx -Wall -Wextra -Wimplicit"

$CC $CFLAGS -o lodepng.o include/lodepng/lodepng.c
$CC $CFLAGS filters.c
$CC $CFLAGS resize.c
$CC $CFLAGS cli.c

$CC -o resize *.o

mkdir -p build
mv resize build

rm *.o
