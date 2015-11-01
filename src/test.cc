#include <stdio.h>
#include "sequences/sequence_file.h"

int main() {
	printf ("[] Compilation of seqlib was successful!\n");
	fflush(stdout);
	printf ("\n");
	fflush(stdout);

	printf ("[] Testing the SequenceFile class.\n");
	fflush(stdout);
	SequenceFile fastqfile("../sample-data/test.fastq");
	fastqfile.Verbose(stdout);
	printf ("\n");

	printf ("[] Done!\n");
	printf ("\n");
	
	return 0;
}
