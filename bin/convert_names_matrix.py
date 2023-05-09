#!/usr/bin/env python3

"""
Change sample names in the STAR counts matrix.
"""

import sys

with open(sys.argv[1], "r") as f:
    x = f.readline()
    z = f.readlines()
    y = []
    for i in x.split("\t"):
        y.append(i)
    f.close()

for k, v in enumerate(y):
    y[k] = v[0:6].replace("-", "_")
s = ""
for i in y:
    s += f"{i}\t"
s = s[0:-1]
s += "\n"

z.pop(0)
z.insert(0, s)

with open(f"mod-{sys.argv[1]}", "w") as f:
    for i in z:
        f.write(i)
    f.close()
