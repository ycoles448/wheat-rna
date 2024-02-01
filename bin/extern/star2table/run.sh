#!/usr/bin/env -S bash

BIN="star2table"
# DEBUG=0

make clean
[ -z "${DEBUG}" ] || echo "Debug enabled"
[ -z "${DEBUG}" ] && make release || make debug

for i in wheat ptr; do
    [ -z "${DEBUG}" ] && \
        "./${BIN}" -o "${i}.tsv" -- ../../../data/star/*${i}*/Log.final.out ||
        valgrind "./${BIN}" -o "${i}.tsv" -- ../../../data/star/*${i}*/Log.final.out
done
