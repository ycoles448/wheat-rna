#ifndef TYPES_H
#define TYPES_H

typedef struct Log {
  int reads_input;
  int mapped_unique;
  int mapped_multi;
  int mapped_excess;
  int unmapped_short;
  int unmapped_mismatch;
  int splices_total;
  int splices_annotated;
  int splices_nc;

  float mapped_avg_length;
} Log;

typedef struct FileList {
  int len;
  char **list;
} FileList;

#endif
