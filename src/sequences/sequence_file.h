/*
 * Copyright 2014, Ivan Sovic.
 * All rights reserved.
 *
 * SequenceFile.h
 *
 *  Created on: 15 May, 2014
 *      Author: Ivan Sovic
 */

#ifndef SEQUENCEFILE_H_
#define SEQUENCEFILE_H_

#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <zlib.h>
//#include "bwa/bwa.h"
//#include "bwa/bwamem.h"
//#include "bwa/kvec.h"
//#include "sequences/utils.h"
#include "sequences/single_sequence.h"
#include "sequences/sequence_alignment.h"
#include "utility/utility_conversion-inl.h"

#include "sequences/kseq.h"

#define SeqFmtToString(x) (x == "fasta" || x == "fastq" || x == "fa" || x == "fq" || x == "fasta.gz" || x == "fastq.gz" || x == "fa.gz" || x == "fq.gz") ? (SEQ_FORMAT_FASTQ) : \
        (x == "gfa" || x == "gfa.gz") ? (SEQ_FORMAT_GFA) : \
        (x == "sam" || x == "sam.gz") ? (SEQ_FORMAT_SAM) : (x == "auto") ? (SEQ_FORMAT_AUTO) : (SEQ_FORMAT_UNKNOWN)

enum SequenceFormat {
  SEQ_FORMAT_UNKNOWN = 0,
  SEQ_FORMAT_AUTO = 1,
  SEQ_FORMAT_FASTQ = 2,     // Reads either FASTA or FASTQ.
  SEQ_FORMAT_GFA = 3,       // Reads the Graphical Assembly Format
  SEQ_FORMAT_SAM = 4,        // Reads the Sequence Alignment/Mapping Format.
//  SEQ_FORMAT_BAM = 5      // TODO: Not implemented yet.
};

//KSEQ_DECLARE(gzFile)
KSEQ_INIT(gzFile, gzread)

// Type defined to hold sequences in the SequenceFile container.
typedef std::vector<SingleSequence *> SequenceVector;

struct seq_sort_key {
  inline bool operator() (const SingleSequence* op1, const SingleSequence* op2) {
    if ((op1)->get_aln().pos != (op2)->get_aln().pos) {
      return ((op1)->get_aln().pos < (op2)->get_aln().pos);
    }
    return (std::string(op1->get_header()) < std::string(op2->get_header()));
  }
};

// This is used for parsing and storing sequences. Supports
// parsing from FASTA and FASTQ files. Files can be parsed
// entirely, or in batches, where the size of a batch can
// be limited either by the number of sequences to load, or
// by the size in MegaBytes to occupy in memory.
// Sample usage 1:
//    SequenceFile sequences1("test1.fasta");
//    sequences1.Verbose(stdout);
//
// Sample usage 2:
//    SequenceFile sequences2;
//    sequences2.OpenFileForBatchLoading("test2.fasta");
//    sequences2.LoadNextBatchNSequences(15);
//    sequences2.LoadNextBatchInMegabytes(2048);
//    sequences2.CloseFileAfterBatchLoading();
//
class SequenceFile {
 public:
  SequenceFile();

  // Initializes the object, and calls LoadAll, with the extension of the file_path as the in format.
  SequenceFile(std::string file_path);

//  // Initializes the object, and calls LoadAll.
//  SequenceFile(SequenceFormat seq_file_fmt, std::string file_path);

  // Initializes the object, and calls OpenFileForBatchLoading and
  // LoadNextBatchNSequences respectively if num_seqs_to_load > 0, or LoadAll otherwise.
  SequenceFile(SequenceFormat seq_file_fmt, std::string file_path,bool convert_to_uppercase=true, uint64_t num_seqs_to_load=0);

  ~SequenceFile();

  // Sets all member variables to initial values, and frees any memory
  // if it was allocated, followed by setting all pointers to NULL value.
  void Clear();

  // Clears only data related to sequences (concretely, the sequences_
  // vector and the current_data_size_ value).
  void ClearOnlyData();

  // Adds a new sequence to the container.
  // Inputs:
  //    sequence  - Pointer to a SingleSequence object holding a new sequence.
  void AddSequence(SingleSequence *sequence, bool needs_destruction);

  // Given the path to a FASTA or a FASTQ file, this function opens the file
  // handlers to enable loading of data from the file.
  // Input:
  //    file_path - Path to the FASTA/FASTQ file.
  // Return:
  //    Returns 0 if successful.
  int OpenFileForBatchLoading(std::string file_path);

  // This function closes the file handles opened with OpenFileForBatchLoading
  // function. Must always be used in pair with OpenFileForBatchLoading.
  // Return:
  //    Returns 0 if successful.
  int CloseFileAfterBatchLoading();

  // Given the path to a file, this function loads all
  // sequences present in the file into this object.
  // Input:
  //    file_path - Path to the FASTA/FASTQ/GFA/... file.
  // Return:
  //    Returns 0 if successful.
  int LoadAll(SequenceFormat seq_file_fmt, std::string file_path, bool convert_to_uppercase=true, bool randomize_non_acgt_bases=false);
  int LoadAllAsBatch(SequenceFormat seq_file_fmt, bool convert_to_uppercase=true, bool randomize_non_acgt_bases=false);

  // This function loads only a part of the sequences present in the
  // file into this object. Before using LoadNextBatchNSequences a file
  // needs to explicitly be opened with the OpenFileForBatchLoading function.
  // This function loads N sequences at once. When called the next time,
  // it will load another batch of N sequences, until EOF is reached.
  // After all batches have been loaded from file, users must call
  // the CloseFileAfterBatchLoading function themselves.
  // Input:
  //    num_seqs_to_load  - Number of sequences to load in a single batch.
  //                        If bigger than total number of sequences in the
  //                        input file, the rest of the file until EOF
  //                        will be loaded.
  //                        If 0, the entire dataset will be loaded.
  // Return:
  //    Returns 0 if successful, -1 if no more sequences can be loaded
  //    (i.e. EOF), and otherwise if unsuccessful.
  int LoadNextBatchNSequences(SequenceFormat seq_file_fmt, uint64_t num_seqs_to_load, bool convert_to_uppercase=true, bool randomize_non_acgt_bases=false);

  // This function loads only a part of the sequences present in the
  // file into this object. Before using LoadNextBatchInMegabytes a file
  // needs to explicitly be opened with the OpenFileForBatchLoading function.
  // This function loads only a fixed amount of MB of sequences at once.
  // When called the next time, it will load another batch fixed size,
  // until EOF is reached.
  // After all batches have been loaded from file, users must call
  // the CloseFileAfterBatchLoading function themselves.
  // Input:
  //    megabytes_to_load - Maximum amount of memory (in MB) to load in a single
  //                        batch. The loaded batch will exceed the value of
  //                        this parameter by the size of the last loaded
  //                        sequence. If the parameter's value is bigger
  //                        than the file size, the entire file will be loaded.
  //                        If 0, the entire dataset will be loaded.
  // Return:
  //    Returns 0 if successful, -1 if no more sequences can be loaded
  //    (i.e. EOF), and otherwise if unsuccessful.
  int LoadNextBatchInMegabytes(SequenceFormat seq_file_fmt, uint64_t megabytes_to_load, bool convert_to_uppercase=true, bool randomize_non_acgt_bases=false);

  // Calculates the size of the sequences that it currently occupies in memory.
  // Calculation of size includes the sum of lengths of the headers, the data
  // and the quality scores of all sequences. It does not include other member
  // variables of this class.
  // Inputs:
  //    memory_unit - specifies the memory unit in which the return value
  //                  should be given (byte, kB, MB or GB). These are
  //                  predefined as constants: MEMORY_UNIT_BYTE,
  //                  MEMORY_UNIT_KILOBYTE, MEMORY_UNIT_MEGABYTE and
  //                  MEMORY_UNIT_GIGABYTE.
  // Return:
  //    Returns the size of the sequences occupying this container.
  uint64_t CalculateTotalSize(int32_t memory_unit=MEMORY_UNIT_BYTE);

  // Returns the sum of number of bases of all the sequences in the file.
  uint64_t GetNumberOfBases();

  // Outputs the contents of this object to the stream given by file pointer.
  // Inputs:
  //    fp  - file pointer to an open file. Can also be stdout and stderr.  void Verbose(FILE *fp);
  void Verbose(FILE *fp);

  // Sorts sequences ascending in this order: 1. position
  void Sort();

  const SequenceVector& get_sequences() const;
  void set_sequences(const SequenceVector& sequences);

  uint64_t get_current_batch_id() const;
  void set_current_batch_id(uint64_t currentBatchId);
  uint64_t get_current_batch_starting_sequence_id() const;
  void set_current_batch_starting_sequence_id(uint64_t currentBatchStartingSequenceId);

 private:
  SequenceVector sequences_;  // Vector holding all the sequences in the file (or in a batch).
  std::vector<bool> destroy_seq_; // Since sequences can be added via AddSequence function from the outside, automatic destruction of them could lead to a segfault or a double free problem. This array keeps track which sequences need to be destroyed by the destructor.
  std::string open_file_path_;  // Path to the sequences file that is currently opened (i.e. during batch loading).
  SequenceFormat seq_file_fmt_;
  kseq_t *bwa_seq_;  // Variable used by BWA's functions for file parsing.
  gzFile gzip_fp_;  // File pointer for an opened FASTA/FASTQ file.
  uint64_t current_batch_id_;  // ID of the current batch (ordinal number).
  uint64_t current_batch_starting_sequence_id_;  // Absolute ID of the first sequence in this object. If all sequences loaded at once is equal to 0, otherwise to the absolute ID of the starting sequence in the batch.
  uint64_t current_data_size_;  // When new sequences are added to this object
                                // using the AddSequence, their size (in bytes)
                                // is automatically calculated. Used for batch
                                // loading of fixed size of sequences.

  int LoadSeqs_(SequenceFormat seq_file_fmt, int64_t num_seqs_to_load, int64_t megabytes_to_load, bool convert_to_uppercase=true, bool randomize_non_acgt_bases=false);

  int LoadSeqsFromFastq_(int64_t num_seqs_to_load, int64_t megabytes_to_load, bool convert_to_uppercase=true, bool randomize_non_acgt_bases=false);

  // Given the path to a GFA file, this function loads all
  // sequences present in the file (lines beginning with 'S') into this object.
  // Input:
  //    file_path - Path to the GFA file.
  // Return:
  //    Returns 0 if successful.
  int LoadSeqsFromGFA_(int64_t num_seqs_to_load, int64_t megabytes_to_load, bool convert_to_uppercase=true, bool randomize_non_acgt_bases=false);

  // Given the path to a SAM file, this function loads all
  // sequences and their alignments present in the file into this object.
  // Input:
  //    file_path - Path to the SAM file.
  // Return:
  //    Returns 0 if successful.
  int LoadSeqsFromSAM_(int64_t num_seqs_to_load, int64_t megabytes_to_load, bool convert_to_uppercase=true, bool randomize_non_acgt_bases=false);

  // Reads a string line from a plain/.gz file. Lines are terminated by '\n' or '\0' characters.
  int ReadGZLine_(gzFile gzip_fp, std::string &ret);

  std::string GetFileExt_(std::string path);
};

#endif /* SEQUENCEFILE_H_ */
