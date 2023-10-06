#ifndef FILES_H
#define FILES_H

#include <stdio.h>

#include "types.h"

void stripSpaceLeft(char *s);
FileList getFileList(int argc, char **argv);
void writeTable(int len, Log *logs);
Log readLog(FILE *f);
void printLog(Log log);

#endif
