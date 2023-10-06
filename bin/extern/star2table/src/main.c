#include <stdio.h>
#include <stdlib.h>

#include "files.h"
#include "macros.h"

int main(int argc, char **argv) {
  FILE *f;

  FileList files;
  Log *logs;

  int i;

  /* Handle arguments */
  if (argc < 1) {
    printf(PROMPT_HELP, argv[0]);
    return ERROR_HELP;
  };

  files = getFileList(argc, argv);
  for (i = 0; i < files.len; i++) {
    printf("%s\n", files.list[i]);
  }

  for (i = 0; i < files.len; i++) {
    /* Load file */
    f = fopen(files.list[i], "r");
    if (f == NULL) {
      perror("main.c:fopen()	" PROMPT_ERROR_FILE);
      return ERROR_FILE;
    }

    /* Read in logs */
    printf("Reading log file:	%s\n", files.list[i]);
    logs = (Log *)malloc(files.len * sizeof(Log));

    logs[i] = readLog(f);
    printLog(logs[i]);

    fclose(f);
  }

  /* Write logs to file */
  writeTable(files.len, logs);

  /* Free memory */
  free(files.list);
  free(logs);

  return 0;
}
