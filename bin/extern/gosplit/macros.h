#ifndef MACROS_H
#define MACROS_H

#define READ_BUFFER 4095

#define PROMPT_HELP \
    "%s go_file [output_prefix] [list_namespace]\n" \
    "\n" \
    " go_file: GO database in OBO format \n" \
    " output_prefix: Prefix for output files, defaults to map-{namespace}.tsv (optional)\n" \
    " list_namespace: A comma-separated list of namespaces within the GO database\n"
#define PROMPT_ERROR_PREFIX "[ERROR]"

#define STRING_ID "id: "
#define STRING_NAME "name: "
#define STRING_NS "namespace: "

#endif
