#!/usr/bin/env python3

"""
Reduce data from experiment metadata to match reads.
"""

import os
import sys
import re

files = []

# Get flags
flags = []
for i, arg in enumerate(sys.argv):
    # Short opts
    if (re.match("^-", arg)) and (arg[0:2] != "--"):
        flag = arg[1:]
        flags.append((i, flag))

    # Long opts
    elif re.match("^--", arg):
        flag = arg[2:]
        flags.append((i, flag))

for flag in flags:
    if flag[1] == "i":
        infile = sys.argv[flag[0] + 1]
        continue
    elif flag[1] == "o":
        outfile = sys.argv[flag[0] + 1]
        continue

# Check flags for duplicates
dups = []
for i, flag in enumerate(flags):
    for j, flag in enumerate(flags):
        if i < j:
            if flags[i] == flags[j]:
                print(f"Found flag duplicates at {i}:{flags[i]} and {j}:{flags[j]}")
                dups.append((j, flags[j]))


print(f"Provided flags: {flags}")
if len(dups) > 0:
    print("Flag duplicates found, exiting...")
    sys.exit(1)

data = []
with open(infile, "r") as f:
    buffer = f.readlines()
    mem = ""

    for i, line in enumerate(buffer[1:]):
        x = line.split("\t")
        unit = x[0]
        row = x[1]
        range = x[2]
        treatment = x[3]
        pool = x[4]
        time = x[5]
        control = x[6]
        cultivar = x[7]
        pot = x[8]
        sample = x[9].strip()
        string = (
            f"{sample}\t{cultivar}\t{control}\t{time}\t{treatment}\t{pool}\t{row}\n"
        )

        if mem != sample:
            data.append(string)
        mem = sample

with open(outfile, "w") as f:
    f.write("sample\tcultivar\tcontrol\ttime\ttreatment\tpool\trow\n")
    for i in data:
        f.write(i)
    f.close()
