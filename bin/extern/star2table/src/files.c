#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "macros.h"
#include "types.h"

void stripSpaceLeft(char *s) {
  char *d = s;

  while (*d == ' ') {
    d++;
  }

  while (*d != '\n') {
    /* printf("%c", *d); */
    *s = *d;
    s++, d++;
  }

  /* Fill end of string with spaces */
  while (*s != '\n') {
    *s = ' ';
    s++;
  }
}

FileList getFileList(int argc, char **argv) {
  FileList files = {0, NULL};

  int i, start;

  /* Get list of files to process */
  /* Check for "--" or lack of flags */
  for (i = 0; i < argc; i++) {
    if (strcmp(argv[i], "--") == 0) break;
  }

  start = i + 1;
  files.len = argc - i - 1;
  files.list = (char **)malloc(files.len * sizeof(char *));

#ifdef DEBUG
  printf("len:	%i\n", files.len);
  printf("start: %i\n", start);
  printf("\n");
#endif
  for (i = 0; i < files.len; i++) {
    files.list[i] = argv[start + i];

#ifdef DEBUG
    printf("i:	%i\n", i);
    printf("file:	%s\n", files.list[i]);
#endif
  }

  return files;
}

void writeTable(int len, Log *logs) {
  int i;

  FILE *f = fopen("out.txt", "w");
  char *headerNames[NSTATS];
  char *header;

  /* Allocate memory */
  for (i = 0; i < NSTATS; i++) {
    headerNames[i] = (char *)malloc(READ_BUFFER_LINE * sizeof(char));
  }
  header = (char *)malloc(READ_BUFFER_LINE * sizeof(char));

  headerNames[0] = TABLE_INPUT_READS;
  headerNames[1] = TABLE_MAPPED_UNIQUE;
  headerNames[2] = TABLE_MAPPED_MULTI;
  headerNames[3] = TABLE_MAPPED_EXCESS;
  headerNames[4] = TABLE_UNMAPPED_SHORT;
  headerNames[5] = TABLE_UNMAPPED_MISMATCH;
  headerNames[6] = TABLE_SPLICES_TOTAL;
  headerNames[7] = TABLE_SPLICES_ANNOTATED;
  headerNames[8] = TABLE_SPLICES_NC;

  /* for (i = 0; i < NSTATS; i++) { */
  /*   header = strcat(header, headerNames[i]); */
  /*   header = strcat(header, "	"); */
  /* } */
  printf("%s\n", header);

  /* Write header */
  /* fputs(); */

  /* Write sample data */
  for (i = 0; i < len; i++) {
    /* logs[i]; */
  }

  /* Free memory */
  for (i = 0; i < 9; i++) {
    free(headerNames[i]);
  }
  free(header);

  fclose(f);
}

Log readLog(FILE *f) {
  Log log = {};
  char *str = (char *)malloc(READ_BUFFER_LINE * sizeof(char));

  while (fgets(str, READ_BUFFER_LINE, f) != NULL) {
    stripSpaceLeft(str);
    sscanf(str, STR_INPUT_READS, &(log.reads_input));
    sscanf(str, STR_MAPPED_UNIQUE, &(log.mapped_unique));
    sscanf(str, STR_MAPPED_MULTI, &(log.mapped_multi));
    sscanf(str, STR_MAPPED_EXCESS, &(log.mapped_excess));
    sscanf(str, STR_UNMAPPED_SHORT, &(log.unmapped_short));
    sscanf(str, STR_UNMAPPED_MISMATCH, &(log.unmapped_mismatch));
    sscanf(str, STR_SPLICES_TOTAL, &(log.unmapped_mismatch));
    sscanf(str, STR_SPLICES_ANNOTATED, &(log.splices_annotated));
    sscanf(str, STR_SPLICES_NC, &(log.splices_nc));
  }

  free(str);

  return log;
}

void printLog(Log log) {
  printf(PROMPT_INPUT_READS, log.reads_input);
  printf(PROMPT_MAPPED_UNIQUE, log.mapped_unique);
  printf(PROMPT_MAPPED_MULTI, log.mapped_multi);
  printf(PROMPT_MAPPED_EXCESS, log.mapped_excess);
  printf(PROMPT_UNMAPPED_SHORT, log.unmapped_short);
  printf(PROMPT_UNMAPPED_MISMATCH, log.unmapped_mismatch);
  printf(PROMPT_SPLICES_TOTAL, log.splices_total);
  printf(PROMPT_SPLICES_ANNOTATED, log.splices_annotated);
  printf(PROMPT_SPLICES_NC, log.splices_nc);
}
