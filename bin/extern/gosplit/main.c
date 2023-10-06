#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "macros.h"

int main(int argc, char **argv) {
    FILE *f;

    char *str;
    char *id;
    char *name;

    int i;
    int ret;

    if (argc < 2) {
        printf(PROMPT_HELP, argv[0]);
        return 1;
    }

    /* Read in database */
    f = fopen(argv[1], "r");

    if (f == NULL) {
        perror(PROMPT_ERROR_PREFIX);
        return -1;
    }

    id = malloc(READ_BUFFER * sizeof(char));
    name = malloc(READ_BUFFER * sizeof(char));

    /* Readlines for the namespace flag */
    str = (char *)malloc(READ_BUFFER * sizeof(char));

    /* Read file in chunks */
    while (fgets(str, READ_BUFFER, f) != NULL) {

        /* Check id */
        if (strstr(str, STRING_ID) != NULL) {
            printf("id: %s", str);
        }

        /* Check name */
        if (strstr(str, STRING_NAME) != NULL) {
            printf("id: %s", str);
        }

        /* Check namespace */
        if (strstr(str, STRING_NS) != NULL) {
            printf("name: %s", str);
        }

        /* printf("%s", str); */
    }
    printf("\n");

    /* Free memory */
    free(id);
    free(name);
    free(str);

    return 0;
}
