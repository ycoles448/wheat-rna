#!/usr/bin/env python3

"""
Extract GO terms using an Interpro tsv output and the reference GFF3 file
containing CDS annotations, taking the path to the tsv file and GFF3 file
as arguments, returning the output as `OUTPUT.map`.
"""

import sys

from pprint import pprint

FILE = sys.argv[1]
REF = sys.argv[2]

# Load Interpro data
with open(FILE, "r") as f:
    buffer = f.readlines()

    for i in range(0, 4):
        buffer.pop(0)

    a = {}
    for k, v in enumerate(buffer):
        x = v.split("\t")[0].strip()
        y = v.split("\t")[-1].strip()

        if y[0:2] == "GO":
            a[x] = {}
            a[x]["go"] = []

    for k, v in enumerate(buffer):
        x = v.split("\t")[0].strip()
        y = v.split("\t")[-1].strip()

        if y[0:2] == "GO":
            y = y.split("|")
            a[x]["go"] = a[x]["go"] + y

    f.close()


# # Load gene IDs from reference
# l = set()
# with open(REF, "r") as f:
#     buffer = f.readlines()

#     # Remove header
#     for i in range(0, 7):
#         buffer.pop(0)

#     for k1, v1 in enumerate(buffer):
#         if v1[0] != "#":
#             t = v1.split("\t")[2]
#             z = v1.split("\t")[8].strip().split(";")

#             if t.lower() == "mrna":
#                 try:
#                     x = z[0][3:]
#                     n = z[1][7:]
#                     l.add((x, n))
#                 except IndexError:
#                     print(f"Unable to determine index for line {k}")
#                     exit(1)

#         # print(f"Processing: {v1}")
#         # for k2, v2 in enumerate(a):
#         #     if v1 == v2:
#         #         print(f"Match for ID: {v1}")
#         #         break

#     f.close()

# c = 0
# for k1 in l:
#     for k2 in a:
#         if k1[0].lower() == k2.lower():
#             a[k1[0]]["name"] = k1[1]
#             c += 1
#             break

pprint(a)
# pprint(l)


# Write map to file
with open("OUTPUT.map", "w") as f:
    for v in a.items():
        print(v)
        s = ""
        for i in v[1]["go"]:
            s += f"{i}, "
        s = f"{v[0][0:-2]}\t{s[0:-2]}\n"
        f.write(s)
    f.close()
