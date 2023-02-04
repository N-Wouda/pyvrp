#!/usr/bin/env bash

# This script builds and installs the pyvrp project in the local poetry 
# environment. Any arguments are passed to the compilation step.

set -e;  # exit immediately if any command fails

BUILD_DIR="build";

# Command line argument for buildtype. Supports a few options, including 
# 'debug' and 'release'. If not provided, a debug build is generated.
BUILD_TYPE=${1:-debug};

if [ ! -d "build" ];  # does a build directory already exist?
then
    poetry run meson setup "$BUILD_DIR" --buildtype "$BUILD_TYPE";
else
    poetry run meson setup "$BUILD_DIR" --reconfigure --buildtype "$BUILD_TYPE";
fi

poetry run meson compile -C "$BUILD_DIR";
poetry run meson install -C "$BUILD_DIR";
