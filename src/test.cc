/// This will only get compiled if -DTEST_SEQ_LIB_ is passed to the compiler.
/// Otherwise, the test file will be skipped, not causing compile time errors for the program including this library.
#ifdef TEST_SEQ_LIB_

#include <stdio.h>
#include "sequences/sequence_file.h"
#include "log_system/log_system.h"
#include "utility/utility_general.h"

int main() {
	printf ("[] Compilation of seqlib was successful!\n");
	fflush(stdout);
	printf ("\n");
	fflush(stdout);

	printf ("[] Testing the SequenceFile class. This test should load and parse a sample FASTQ file.\n");
	fflush(stdout);
	SequenceFile fastqfile("sample-data/test.fastq");
	fastqfile.Verbose(stdout);
	printf ("\n");

	printf ("[] Done!\n");
	printf ("\n");

	LogSystem::GetInstance().SetProgramVerboseLevelFromInt(VERBOSE_LEVEL_LOW | VERBOSE_LEVEL_DEBUG);
	LogSystem::GetInstance().WriteLog("Example of how to directly write to a log file!", true);
	LogSystem::GetInstance().WriteLog("This gets written to file but not to stderr, but only in case (LogSystem::LOG_VERBOSE_TYPE & LOG_VERBOSE_STD) == 0.", false);
	LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Example of how to write to stderr."), "LogExample");
	LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED, true, FormatString("This message will only get displayed on a medium verbose level.\n"), "LogExample");
	LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH, true, FormatString("This message will only get displayed on a high verbose level.\n"), "LogExample");
	LogSystem::GetInstance().Log(VERBOSE_LEVEL_LOW | VERBOSE_LEVEL_DEBUG, true, FormatString("This message will only get displayed on a low debug verbose level.\n"), "LogExample");
	LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED | VERBOSE_LEVEL_DEBUG, true, FormatString("This message will only get displayed on a medium debug verbose level.\n"), "LogExample");
	LogSystem::GetInstance().Error(SEVERITY_INT_INFO, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "This is an example of an info message."));
	LogSystem::GetInstance().Error(SEVERITY_INT_WARNING, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "This is an example of a warning message."));
	/// This is commented out, because it would exit the program.
	// LogSystem::GetInstance().Error(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "This is an example of a warning message."));
	// LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "This is an example of a memory exception."));

	SequenceFile samfile("sample-data/test.sam");
	samfile.Verbose(stdout);
	for (int64_t i=0; i<samfile.get_sequences().size(); i++) {
		printf ("[%ld] Pos: %d, Qname: '%s'\n", i, samfile.get_sequences()[i]->get_aln().get_pos(), samfile.get_sequences()[i]->get_header());
	}
	printf ("\n");

	printf ("Sorting the sequences.\n");
	samfile.Sort();
	for (int64_t i=0; i<samfile.get_sequences().size(); i++) {
		printf ("[%ld] Pos: %d, Qname: '%s'\n", i, samfile.get_sequences()[i]->get_aln().get_pos(), samfile.get_sequences()[i]->get_header());
	}
	printf ("\n");
	// samfile.Sort();

	printf ("Testing the correctnes of leading/trailing clipped bases:\n");
	for (int64_t i=0; i<samfile.get_sequences().size(); i++) {
		std::string cigar = samfile.get_sequences()[i]->get_aln().GetCigarString();
		printf ("[%ld] front: %ld, back: %ld, CIGAR start: %s, CIGAR end: %s\n", i, samfile.get_sequences()[i]->get_aln().GetClippedBasesFront(), samfile.get_sequences()[i]->get_aln().GetClippedBasesBack(), GetSubstring((char *) cigar.c_str(), 10).c_str(), GetSubstring((char *) (cigar.c_str() + cigar.size() - 10), 10).c_str() );
	}

	return 0;
}

#endif
