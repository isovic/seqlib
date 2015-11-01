#include <stdio.h>
#include "sequences/sequence_file.h"
#include "log_system/log_system.h"

int main() {
	printf ("[] Compilation of seqlib was successful!\n");
	fflush(stdout);
	printf ("\n");
	fflush(stdout);

	printf ("[] Testing the SequenceFile class. This test should load and parse a sample FASTQ file.\n");
	fflush(stdout);
	SequenceFile fastqfile("../sample-data/test.fastq");
	fastqfile.Verbose(stdout);
	printf ("\n");

	printf ("[] Done!\n");
	printf ("\n");

	LogSystem::GetInstance().SetProgramVerboseLevelFromInt(VERBOSE_LEVEL_LOW | VERBOSE_LEVEL_DEBUG);
	LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Example of how to write to a log file!\n"), "LogExample");
	LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED, true, FormatString("This message will only get displayed on a medium verbose level.\n"), "LogExample");
	LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH, true, FormatString("This message will only get displayed on a high verbose level.\n"), "LogExample");
	LogSystem::GetInstance().Log(VERBOSE_LEVEL_LOW | VERBOSE_LEVEL_DEBUG, true, FormatString("This message will only get displayed on a low debug verbose level.\n"), "LogExample");
	LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED | VERBOSE_LEVEL_DEBUG, true, FormatString("This message will only get displayed on a medium debug verbose level.\n"), "LogExample");
	LogSystem::GetInstance().Error(SEVERITY_INT_INFO, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "This is an example of an info message."));
	LogSystem::GetInstance().Error(SEVERITY_INT_WARNING, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "This is an example of a warning message."));
	LogSystem::GetInstance().Error(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "This is an example of a warning message."));
	/// This is commented out, because it would exit the program.
	/// LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "This is an example of a memory exception."));

	return 0;
}
