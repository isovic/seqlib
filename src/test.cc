#include <stdio.h>
#include "sequences/sequence_file.h"
#include "log_system/log_system.h"

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

	LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "This is an example of a memory exception."));

	return 0;
}
