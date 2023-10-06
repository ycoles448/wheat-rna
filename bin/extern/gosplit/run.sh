#!/usr/bin/env bash

BUILD="debug"
BIN="gosplit"

make clean
[[ "${BUILD}" == "debug" ]] && make debug || make release

./"${BIN}" "${@}"
