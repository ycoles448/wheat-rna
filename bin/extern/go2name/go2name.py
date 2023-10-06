#!/usr/bin/env python3

import sys
import re
import json

map = "map.txt"

def showHelp() -> None:
    print("""
./go2name map_file go_file ontology [annot_file]

    map_file is a file containing tab separated columns of a singular gene ID, and a comma-separated list of GO terms
    go_file is a GO annotation file, in OBO format
    ontology is the namespace to export from the GO database
    annot_file is the output, defaults to "map.txt"
    """)

def __main__() -> int:
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        showHelp()
        return 1

    mapFile = open(sys.argv[1], "r")
    goFile = open(sys.argv[2], "r")

    ontology = sys.argv[3]
    buffer = mapFile.readlines()
    mapFile.close()

    a = {}
    for v in buffer:
        a[v.split("\t")[0].strip()] = dict(id = [i.strip() for i in v.split("\t")[1].split(",")], names = [])

    buffer = goFile.readlines()
    goFile.close()

    b = {}
    for k, v in enumerate(buffer):
        if re.match("id: GO:.*", v):
            i = v.split(" ")[1].strip()
            b[i] = dict(name = buffer[k + 1].split(": ")[1].strip())
            b[i]["ontology"] = buffer[k + 2].split(": ")[1].strip()
            b[i]["genes"] = []

    # Match GO terms to an annotation name
    for v1 in a:
        for i in a[v1]["id"]:
            for v2 in b:
                if i == v2 and b[v2]["ontology"] == ontology:
                    a[v1]["names"].append(b[v2]["name"])
                    break
        # print(f"{v1}: {a[v1]}")

    # Export new map
    if len(sys.argv) == 5:
        outFile = open(sys.argv[4] + f"-{ontology}", "w")
    else:
        outFile = open(f"map-{ontology}.txt", "w")

    str = []
    for v in a:
        strNames = ""
        for i in a[v]["names"]:
            strNames += f"{i};"
        strNames = strNames[0:len(strNames) - 1]
        if len(strNames) == 0:
            continue

        # strIDs = ""
        # for i in a[v]["id"]:
        #     strIDs += f"{i}, "
        # strIDs = strIDs[0:len(strIDs) - 2]

        str.append(f"{v}\t{strNames}\n")
    outFile.writelines(str)

    outFile.close()

    return 0

__main__()
