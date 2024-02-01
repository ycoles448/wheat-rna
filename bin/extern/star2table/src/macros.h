/* Local Variables: */
/* End: */

#ifndef MACROS_H
#define MACROS_H

/* Global variables */
extern char *outfile;

/* Number of individual statistics used from each line */
#define NSTATS 9

#define PROMPT_HELP "Usage: %s [--] files\n"
#define PROMPT_INPUT_READS "Input reads: %i\n"
#define PROMPT_MAPPED_UNIQUE "Mapped reads (unique): %i\n"
#define PROMPT_MAPPED_MULTI "Mapped reads (multi): %i\n"
#define PROMPT_MAPPED_EXCESS "Mapped reads (excess): %i\n"
#define PROMPT_UNMAPPED_SHORT "Unmapped reads (too short): %i\n"
#define PROMPT_UNMAPPED_MISMATCH "Unmapped reads (mismatch): %i\n"
#define PROMPT_SPLICES_TOTAL "Splices (total): %i\n"
#define PROMPT_SPLICES_ANNOTATED "Splices (annotated): %i\n"
#define PROMPT_SPLICES_NC "Splices (non-canonical): %i\n"

#define STR_SEP "	"

#define STR_INPUT_READS "Number of input reads |	%i"
#define STR_MAPPED_UNIQUE "Uniquely mapped reads number |	%i"
#define STR_MAPPED_MULTI "Number of reads mapped to multiple loci |	%i"
#define STR_MAPPED_EXCESS "Number of reads mapped to too many loci |	%i"
#define STR_UNMAPPED_SHORT "Number of reads unmapped: too short |	%i"
#define STR_UNMAPPED_MISMATCH \
  "Number of reads unmapped: too many mismatches |	%i"
#define STR_SPLICES_TOTAL "Number of splices: Total |	%i"
#define STR_SPLICES_ANNOTATED "Number of splices: Annotated (sjdb) |	%i"
#define STR_SPLICES_NC "Number of splices: Non-canonical |	%i"

#define TABLE_HEADER "sample"
#define TABLE_INPUT_READS "input_reads"
#define TABLE_MAPPED_UNIQUE "mapped_unique"
#define TABLE_MAPPED_MULTI "mapped_multiple"
#define TABLE_MAPPED_EXCESS "mapped_multiple_excess"
#define TABLE_UNMAPPED_SHORT "unmapped_short"
#define TABLE_UNMAPPED_MISMATCH "unmapped_mismatch"
#define TABLE_SPLICES_TOTAL "splices_total"
#define TABLE_SPLICES_ANNOTATED "splices_annotated"
#define TABLE_SPLICES_NC "splices_noncanonical"

#define ERROR_FILE 10
#define PROMPT_ERROR_FILE "Could not read input file"

#define ERROR_HELP 1

#define LOG_FILE "Log.out.final"

/* Number of characters per line */
#define READ_BUFFER_LINE 1024

#endif
