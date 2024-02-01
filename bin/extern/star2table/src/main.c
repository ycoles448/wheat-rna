#include <stdio.h>
#include <stdlib.h>

#include "files.h"
#include "macros.h"

/* Initialise global variables */
char *outfile = NULL;

int main(int argc, char **argv) {
  FILE *f;

  FileList files;
  Log *logs;

  int i;

  /* Handle arguments */
  processArgs(argc, argv);

  if (argc < 1) {
    printf(PROMPT_HELP, argv[0]);
    return ERROR_HELP;
  };

  files = getFileList(argc, argv);

  logs = malloc(files.len * sizeof(Log));
  for (i = 0; i < files.len; i++) {
    /* Load file */
    f = fopen(files.list[i], "r");
    if (f == NULL) {
      perror("main.c:fopen()	" PROMPT_ERROR_FILE);
      return ERROR_FILE;
    }

    /* Read in logs */
    /* printf("Reading log file:	%s\n", files.list[i]); */
    logs[i] = readLog(f);
    /* printLog(logs[i]); */

    fclose(f);
  }

  /* Write logs to file */
  writeTable(files.len, logs);

  /* Free memory */
  free(files.list);
  free(logs);

  return 0;
}
