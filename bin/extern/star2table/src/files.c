#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "macros.h"
#include "types.h"

/* Process command line arguments */
void processArgs(int argc, char **argv) {
  int i;
  char *s;

  for (i = 0; i < argc; i++) {
    s = argv[i];

    if (strcmp(s, "-o") == 0) {
      outfile = argv[i + 1];
    };
  }
}

/* Remove spaces from the LHS of a string */
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
  files.list = malloc(files.len * sizeof(char *));

  for (i = 0; i < files.len; i++) files.list[i] = argv[start + i];

  return files;
}

void writeTable(int len, Log *logs) {
  int i;

  FILE *f = fopen(outfile, "w");
  char *headerNames[NSTATS];
  char *header;
  char *line;
  char *tmp;
  /* char *num; */
  char *sep = "	";

  printf("Writing to file:	%s\n", outfile);

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
  /*   printf("%s\n", headerNames[i]); */
  /* } */

  header = calloc(READ_BUFFER_LINE, sizeof(char));
  strcat(header, TABLE_HEADER);
  strcat(header, sep);
  for (i = 0; i < NSTATS; i++) {
    strcat(header, headerNames[i]);
    strcat(header, sep);
  }
  for (i = strlen(header) - strlen(sep); i < strlen(header); i++) header[i] = 0;

  fputs(header, f);
  fputc('\n', f);

  /* Print lines from each log */
  /* Format string for number of reads */
  tmp = calloc(READ_BUFFER_LINE, sizeof(char));
  strcat(tmp, "sample%i");
  strcat(tmp, sep);

  for (i = 0; i < NSTATS; i++) {
    strcat(tmp, "%i");
    strcat(tmp, sep);
  }
  for (i = strlen(tmp) - strlen(sep); i < strlen(tmp); i++) tmp[i] = 0;

  for (i = 0; i < len; i++) {
    line = calloc(READ_BUFFER_LINE, sizeof(char));

    sprintf(line, tmp, i + 1, logs[i].reads_input, logs[i].mapped_unique,
            logs[i].mapped_multi, logs[i].mapped_excess, logs[i].unmapped_short,
            logs[i].unmapped_mismatch, logs[i].splices_total,
            logs[i].splices_annotated, logs[i].splices_nc);

    /* printf("%s\n", line); */
    fputs(line, f);
    fputc('\n', f);

    free(line);
  }

  /* Free memory */
  free(header);
  free(tmp);
  fclose(f);
}

Log readLog(FILE *f) {
  Log log = {};
  char *str = malloc(READ_BUFFER_LINE * sizeof(char));

  while (fgets(str, READ_BUFFER_LINE, f) != NULL) {
    stripSpaceLeft(str);
    sscanf(str, STR_INPUT_READS, &(log.reads_input));
    sscanf(str, STR_MAPPED_UNIQUE, &(log.mapped_unique));
    sscanf(str, STR_MAPPED_MULTI, &(log.mapped_multi));
    sscanf(str, STR_MAPPED_EXCESS, &(log.mapped_excess));
    sscanf(str, STR_UNMAPPED_SHORT, &(log.unmapped_short));
    sscanf(str, STR_UNMAPPED_MISMATCH, &(log.unmapped_mismatch));
    sscanf(str, STR_SPLICES_TOTAL, &(log.splices_total));
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
