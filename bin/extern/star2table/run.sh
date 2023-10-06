#!/usr/bin/env -S bash

BIN="star2table"
# DEBUG=1

make clean
[ -z "${DEBUG}" ] && make release || make debug
[ -z "${DEBUG}" ] && "./${BIN}" "${@}" || valgrind -s --track-origins=yes --leak-check=full "./${BIN}" "${@}"
