/*
 * Copyright 2014, Ivan Sovic.
 * All rights reserved.
 *
 * SequenceFile.cc
 *
 *  Created on: 15 May, 2014
 *      Author: Ivan Sovic
 */

#include "sequences/sequence_file.h"
#include "log_system/log_system.h"
#include <map>
#include <sstream>
#include <algorithm>

//#include "bwa/kseq.h"
//KSEQ_DECLARE(gzFile)



SequenceFile::SequenceFile() {
  bwa_seq_ = NULL;
  Clear();
}

SequenceFile::SequenceFile(std::string file_path) {
  bwa_seq_ = NULL;
  Clear();

  seq_file_fmt_ = SEQ_FORMAT_AUTO;
  LoadAll(seq_file_fmt_, file_path, true);
}

//SequenceFile::SequenceFile(SequenceFormat seq_file_fmt, std::string file_path) {
//  bwa_seq_ = NULL;
//  Clear();
//  seq_file_fmt_ = seq_file_fmt;
//  LoadAll(seq_file_fmt, file_path, true);
//}

SequenceFile::SequenceFile(SequenceFormat seq_file_fmt, std::string file_path, bool convert_to_uppercase, uint64_t num_seqs_to_load) {
  bwa_seq_ = NULL;
  Clear();
  seq_file_fmt_ = seq_file_fmt;
  if (num_seqs_to_load <= 0) {
      LoadAll(seq_file_fmt, file_path, convert_to_uppercase);
  } else {
    OpenFileForBatchLoading(file_path);
    LoadNextBatchNSequences(seq_file_fmt, num_seqs_to_load, convert_to_uppercase);
  }
}

SequenceFile::~SequenceFile() {
  Clear();
}

void SequenceFile::Clear() {
  ClearOnlyData();

  current_batch_id_ = 0;
  current_batch_starting_sequence_id_ = 0;
  open_file_path_ = "";

  if (bwa_seq_ != NULL) {
    kseq_destroy(bwa_seq_);
    bwa_seq_ = NULL;
  }
}

void SequenceFile::ClearOnlyData() {
  for (int64_t i=0; i<sequences_.size(); i++) {
    if (destroy_seq_[i] == true && sequences_[i]) {
      delete sequences_[i];
    }
  }

  sequences_.clear();
  destroy_seq_.clear();
  current_data_size_ = 0;
}

void SequenceFile::AddSequence(SingleSequence *sequence, bool needs_destruction) {
  sequences_.push_back(sequence);
  destroy_seq_.push_back(needs_destruction);
  current_data_size_ += sequence->CalculateTotalSize(kMemoryUnitByte);
}

const SequenceVector& SequenceFile::get_sequences() const {
  return sequences_;
}

uint64_t SequenceFile::GetNumberOfBases() {
  uint64_t ret = 0;
  for (uint64_t i=0; i<sequences_.size(); i++) {
    ret += sequences_[i]->get_sequence_length();
  }
  return ret;
}

void SequenceFile::set_sequences(const SequenceVector& sequences) {
  sequences_ = sequences;
}

int SequenceFile::LoadAll(SequenceFormat seq_file_fmt, std::string file_path, bool convert_to_uppercase, bool randomize_non_acgt_bases) {
  Clear();

  if (OpenFileForBatchLoading(file_path))
    return 1;

  int ret_val = 0;
  if ((ret_val = LoadAllAsBatch(seq_file_fmt, convert_to_uppercase, randomize_non_acgt_bases)) != 0) {
    CloseFileAfterBatchLoading();
    return ret_val;
  }

  return CloseFileAfterBatchLoading();
}

int SequenceFile::OpenFileForBatchLoading(std::string file_path) {
  gzip_fp_ = gzopen(file_path.c_str(), "r");

  if (gzip_fp_ == NULL) {
    FATAL_REPORT(ERR_OPENING_FILE, "File path: '%s'.", file_path.c_str());
    return 1;
  }

  bwa_seq_ = kseq_init(gzip_fp_);

  open_file_path_ = file_path;

  return 0;
}

int SequenceFile::CloseFileAfterBatchLoading() {
  if (bwa_seq_ == NULL) {
    LogSystem::GetInstance().Error(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_CLOSING_FILE, "Offending variable: bwt_seq_."));
    return 1;
  }

  kseq_destroy(bwa_seq_);

  bwa_seq_ = NULL;
  open_file_path_ = "";

  gzclose(gzip_fp_);

  return 0;
}

int SequenceFile::LoadAllAsBatch(SequenceFormat seq_file_fmt, bool convert_to_uppercase, bool randomize_non_acgt_bases) {
  ClearOnlyData();

  current_batch_id_ = 0;
  current_batch_starting_sequence_id_ = 0;

  if (gzip_fp_ == NULL || bwa_seq_ == NULL)
  {
    if (gzip_fp_ == NULL)
      LogSystem::GetInstance().Error(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_OPENING_FILE, "Offending variable: bwt_fp_."));
    if (bwa_seq_ == NULL)
      LogSystem::GetInstance().Error(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_WRONG_PARAMS, "Offending variable: bwt_seq_."));
    return 1;
  }

  int64_t num_loaded_seqs = LoadSeqs_(seq_file_fmt, 0, 0, convert_to_uppercase, randomize_non_acgt_bases);

  // If there are quality values, but for some reason the length of
  // the quality string is different from the length of the sequences
  // then the file is corrupted. Reporting the error here.
  // This check is repeated here because of BWA's kseq_read function
  // which can break on such sequences, thus to report the error non the less.
  // Additionally, if the functionality of kseq_read is changed, the error will
  // be reported anyway.
  if ((bwa_seq_->qual.l > 0 && bwa_seq_->seq.l != bwa_seq_->qual.l) || num_loaded_seqs < 0) {
    ERROR_REPORT(ERR_FILE_DEFORMED_FORMAT, "Quality length not equal to sequence length in FASTQ file '%s'. Sequence ID: %ld. Skipping rest of the file.", open_file_path_.c_str(), sequences_.size());
    return -2;
  }

  // This happens when EOF is reached. Not technically an error, but not as if the batch has been loaded.
  if (sequences_.size() == 0)
    return -1;

  return 0;
}

int SequenceFile::LoadNextBatchNSequences(SequenceFormat seq_file_fmt, uint64_t num_seqs_to_load, bool convert_to_uppercase, bool randomize_non_acgt_bases) {
  current_batch_starting_sequence_id_ += sequences_.size();

  ClearOnlyData();

  if (gzip_fp_ == NULL || bwa_seq_ == NULL)
  {
    if (gzip_fp_ == NULL)
      LogSystem::GetInstance().Error(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_OPENING_FILE, "Offending variable: bwt_fp_."));
    if (bwa_seq_ == NULL)
      LogSystem::GetInstance().Error(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_WRONG_PARAMS, "Offending variable: bwt_seq_."));
    return 1;
  }

  int64_t num_loaded_seqs = 0;
  num_loaded_seqs = LoadSeqs_(seq_file_fmt, num_seqs_to_load, 0, convert_to_uppercase, randomize_non_acgt_bases);

  if ((bwa_seq_->qual.l > 0 && bwa_seq_->seq.l != bwa_seq_->qual.l) || num_loaded_seqs < 0) {
    ERROR_REPORT(ERR_FILE_DEFORMED_FORMAT, "Quality length not equal to sequence length in FASTQ file! Batch ID: %ld. Absolute sequence ID: %ld. Relative sequence ID: %ld. Skipping rest of the file.", current_batch_id_, (current_batch_starting_sequence_id_ + sequences_.size()), sequences_.size());
    CloseFileAfterBatchLoading();
    return -2;
  }

  // This happens when EOF is reached. Not technically an error, but not as if the batch has been loaded.
  if (sequences_.size() == 0)
    return -1;

  return 0;
}

int SequenceFile::LoadNextBatchInMegabytes(SequenceFormat seq_file_fmt, uint64_t megabytes_to_load, bool convert_to_uppercase, bool randomize_non_acgt_bases) {
  current_batch_starting_sequence_id_ += sequences_.size();

  ClearOnlyData();

  if (gzip_fp_ == NULL || bwa_seq_ == NULL)
  {
    if (gzip_fp_ == NULL)
      LogSystem::GetInstance().Error(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_OPENING_FILE, "Offending variable: bwt_fp_."));
    if (bwa_seq_ == NULL)
      LogSystem::GetInstance().Error(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_WRONG_PARAMS, "Offending variable: bwt_seq_."));
    return 1;
  }

  int64_t num_loaded_seqs = 0;
  num_loaded_seqs = LoadSeqs_(seq_file_fmt, 0, megabytes_to_load, convert_to_uppercase, randomize_non_acgt_bases);

  if ((bwa_seq_->qual.l > 0 && bwa_seq_->seq.l != bwa_seq_->qual.l) || num_loaded_seqs < 0) {
    LogSystem::GetInstance().Error(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_DEFORMED_FORMAT, "Quality length not equal to sequence length in FASTQ file! Batch ID: %ld. Absolute sequence ID: %ld. Relative sequence ID: %ld. Skipping rest of the file.", current_batch_id_, (current_batch_starting_sequence_id_ + num_loaded_seqs), num_loaded_seqs));
    CloseFileAfterBatchLoading();
    return -2;
  }

  // This happens when EOF is reached. Not technically an error, but not as if the batch has been loaded.
  if (sequences_.size() == 0)
    return -1;

  return 0;
}

uint64_t SequenceFile::CalculateTotalSize(int32_t memory_unit) {
  uint64_t total_size = 0;

  for (SequenceVector::iterator sequence_iterator = sequences_.begin();
      sequence_iterator != sequences_.end(); sequence_iterator++) {
    total_size += (*sequence_iterator)->CalculateTotalSize(kMemoryUnitByte);
  }

  if (memory_unit == MEMORY_UNIT_BYTE)
    total_size = total_size / ((uint64_t) 1);
  else if (memory_unit == MEMORY_UNIT_KILOBYTE)
    total_size = total_size / ((uint64_t) 1024);
  else if (memory_unit == MEMORY_UNIT_MEGABYTE)
    total_size = total_size / (((uint64_t) 1024) * ((uint64_t) 1024));
  else if (memory_unit == MEMORY_UNIT_GIGABYTE)
    total_size = total_size / (((uint64_t) 1024) * ((uint64_t) 1024) * ((uint64_t) 1024));
  else
    LogSystem::GetInstance().Error(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_WRONG_PARAMS, "Memory unit not recognized! Returning value in bytes."));

  return total_size;
}

void SequenceFile::Verbose(FILE *fp) {
  fprintf (fp, "Num sequences: %ld\n", sequences_.size());
  fprintf (fp, "Currently open file (only when batch loading): '%s'\n", open_file_path_.c_str());
  fprintf (fp, "Current batch ID: %ld\n", current_batch_id_);
  fprintf (fp, "Current batch starting sequence ID: %ld\n", current_batch_starting_sequence_id_);

  fprintf (fp, "\n");

//  for (SequenceVector::iterator sequence_iterator = sequences_.begin();
//      sequence_iterator != sequences_.end(); sequence_iterator++) {
//    (*sequence_iterator)->Verbose(fp);
//    fprintf (fp, "\n");
//  }

  fflush(fp);
}

uint64_t SequenceFile::get_current_batch_id() const {
  return current_batch_id_;
}

void SequenceFile::set_current_batch_id(uint64_t currentBatchId) {
  current_batch_id_ = currentBatchId;
}

uint64_t SequenceFile::get_current_batch_starting_sequence_id() const {
  return current_batch_starting_sequence_id_;
}

void SequenceFile::set_current_batch_starting_sequence_id(uint64_t currentBatchStartingSequenceId) {
  current_batch_starting_sequence_id_ = currentBatchStartingSequenceId;
}

std::string SequenceFile::GetFileExt_(std::string path) {
  int32_t pos = path.find_last_of(".");
  std::string ext = path.substr(pos + 1);
  if (ext == "gz") {
    ext = path.substr(path.find_last_of(".", (pos-1)) + 1);
//    printf ("ext = %s\n", ext.c_str());
  }
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
  return ext;
}

int SequenceFile::LoadSeqs_(SequenceFormat seq_file_fmt, int64_t num_seqs_to_load, int64_t megabytes_to_load, bool convert_to_uppercase, bool randomize_non_acgt_bases) {
  if (seq_file_fmt == SEQ_FORMAT_AUTO) {
    seq_file_fmt = SeqFmtToString(GetFileExt_(open_file_path_));
  }

  if (seq_file_fmt == SEQ_FORMAT_FASTQ) {
    return LoadSeqsFromFastq_(num_seqs_to_load, megabytes_to_load, convert_to_uppercase, randomize_non_acgt_bases);

  } else if (seq_file_fmt == SEQ_FORMAT_GFA) {
    return LoadSeqsFromGFA_(num_seqs_to_load, megabytes_to_load, convert_to_uppercase, randomize_non_acgt_bases);

  } else if (seq_file_fmt == SEQ_FORMAT_SAM) {
    return LoadSeqsFromSAM_(num_seqs_to_load, megabytes_to_load, convert_to_uppercase, randomize_non_acgt_bases);

  } else {
    FATAL_REPORT(ERR_UNEXPECTED_VALUE, "Input sequence file format unknown!\n");
    return -1;
  }

  return 0;
}

int SequenceFile::LoadSeqsFromFastq_(int64_t num_seqs_to_load, int64_t megabytes_to_load, bool convert_to_uppercase, bool randomize_non_acgt_bases) {
  int32_t l;
  uint64_t id = 0;
  uint64_t id_absolute = current_batch_starting_sequence_id_;
  SingleSequence *sequence = NULL;

  while ((l = kseq_read(bwa_seq_)) >= 0) {
    sequence = new SingleSequence();

    // BWA's parsing functions split the headers in two parts, name
    // and comment. We join them here again, but we must check that
    // these are not null pointers or empty strings.
    std::string header("");
    if (bwa_seq_->name.l > 0)
      header += std::string(bwa_seq_->name.s);
    if (bwa_seq_->comment.l > 0)
      header += std::string(" ") + std::string(bwa_seq_->comment.s);

    if (!bwa_seq_->qual.l) {  // If bwa_seq_->qual.l is equal to 0, then we are loading a FASTA file. In this case, only header and data need to be initialized.
      sequence->InitHeaderAndDataFromAscii((char *) header.c_str(),
                                                header.length(),
                                                (int8_t *) bwa_seq_->seq.s,
                                                bwa_seq_->seq.l, id, id_absolute);
    } else {  // If we got here, we are loading a FASTQ file. All three components (header, data and quality scores) need to be initialized.
      if (bwa_seq_->seq.l != bwa_seq_->qual.l) {
        LogSystem::GetInstance().Error(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_DEFORMED_FORMAT, "Quality length not equal to sequence length in FASTQ file! Batch ID: %ld. Absolute sequence ID: %ld. Relative sequence ID: %ld. Skipping rest of the file.", current_batch_id_, (current_batch_starting_sequence_id_ + id), id));
        CloseFileAfterBatchLoading();
        return 1;
      }

      sequence->InitAllFromAscii((char *) header.c_str(), header.length(),
                                      (int8_t *) bwa_seq_->seq.s,
                                      (int8_t *) bwa_seq_->qual.s, bwa_seq_->seq.l, id, id_absolute);
    }

    if (randomize_non_acgt_bases == true) {
      sequence->RandomizeNonACGTBases();
    }
    if (convert_to_uppercase == true) {
      sequence->BasesToUppercase();
    }

    AddSequence(sequence, true);
    id += 1;  // Increment the relative sequence id counter.
    id_absolute += 1;

    if (num_seqs_to_load > 0 && id >= num_seqs_to_load)  // Batch loading stopping condition.
      break;
    if (megabytes_to_load > 0 && ConvertFromBytes(MEMORY_UNIT_MEGABYTE, current_data_size_) >= megabytes_to_load)  // Batch loading stopping condition.
      break;
  }

  if (l == -2)
    return -2;

  return id;
}

int SequenceFile::ReadGZLine_(gzFile gzip_fp, std::string &ret) {
  const int32_t BUFF_SIZE = 4096;
  int8_t buff[BUFF_SIZE + 1];
  std::stringstream ss;
  int32_t ret_len = 0;
  int32_t read_bytes = 0;
  bool ln_end_found = false;
  int32_t ln_end_pos = 0;
  bool is_eof = false;
  while ((read_bytes = gzread(gzip_fp, &buff[0], BUFF_SIZE)) > 0) {
    buff[read_bytes] = '\0';
    for (int32_t i=0; i<read_bytes; i++) {
      if (buff[i] == '\n') {
        buff[i] = '\0';
        gzseek(gzip_fp, -(read_bytes - i - 1), SEEK_CUR);
        ln_end_found = true;
        ln_end_pos = i;
        break;
      }
    }
    std::string buff_str((char *) &buff[0]);
    ss << buff_str;
    ret_len += buff_str.length();

    if (ln_end_found) break;
  }

  ret = ss.str();
  // This allows for empty lines to be output as well (but not if there is an EOF, i.e. read_bytes == 0).
  if (ret_len > 0 || read_bytes > 0) return 0;

  return 1;
}

int SequenceFile::LoadSeqsFromGFA_(int64_t num_seqs_to_load, int64_t megabytes_to_load, bool convert_to_uppercase, bool randomize_non_acgt_bases) {
  uint64_t id = 0;
  uint64_t id_absolute = current_batch_starting_sequence_id_;
  SingleSequence *sequence = NULL;
  std::string line;

  std::map<std::string, SingleSequence *> seq_header_map;

  while (!ReadGZLine_(gzip_fp_, line)) {
    sequence = new SingleSequence();

    std::istringstream ss(line);
    std::string keyword;
    ss >> keyword;

    if (keyword == "S") {
      std::string header;
      std::string seq;
      std::string ln;
      ss >> header >> seq >> ln;

      sequence->InitHeaderAndDataFromAscii((char *) header.c_str(),
                                                header.length(),
                                                (int8_t *) seq.c_str(), seq.length(), id, id_absolute);

      if (randomize_non_acgt_bases == true) {
        sequence->RandomizeNonACGTBases();
      }
      if (convert_to_uppercase == true) {
        sequence->BasesToUppercase();
      }

      AddSequence(sequence, true);
      seq_header_map[header] = sequence;  // Hash the sequence by it's name for faster search when adding other GFA lines.
      id += 1;  // Increment the relative sequence id counter.
      id_absolute += 1;

      if (num_seqs_to_load > 0 && id >= num_seqs_to_load)  // Batch loading stopping condition.
        break;
      if (megabytes_to_load > 0 && ConvertFromBytes(MEMORY_UNIT_MEGABYTE, current_data_size_) >= megabytes_to_load)  // Batch loading stopping condition.
        break;

    } else if (keyword == "a") {
      GFAGoldenPath gp;
      SequenceGFA::ParseGoldenPath(line, gp);

      auto it = seq_header_map.find(gp.utg_name);
      if (it == seq_header_map.end()) {
        FATAL_REPORT(ERR_UNEXPECTED_VALUE, "Unitig name not found in previously loaded sequences! utg_name = '%s'", gp.utg_name.c_str());
      }
      it->second->gfa().AddGoldenPath(gp);

    } else if (keyword == "L") {
      GFAOverlap ov;
      SequenceGFA::ParseOverlap(line, ov);

      auto it1 = seq_header_map.find(ov.seg1_name);
      if (it1 != seq_header_map.end()) {
        it1->second->gfa().AddOverlap(ov);
      }

      auto it2 = seq_header_map.find(ov.seg2_name);
      if (it2 != seq_header_map.end()) {
        it2->second->gfa().AddOverlap(ov);
      }

      if (it1 == seq_header_map.end() || it2 == seq_header_map.end()) {
        FATAL_REPORT(ERR_UNEXPECTED_VALUE, "Overlap cannot be added to GFA info because either seg1_name or seg2_name cannot be found in the sequences loaded so far. line = '%s'", line.c_str());
      }
    }
  }

  return id;
}

int SequenceFile::LoadSeqsFromSAM_(int64_t num_seqs_to_load, int64_t megabytes_to_load, bool convert_to_uppercase, bool randomize_non_acgt_bases) {
  uint64_t id = 0;
  uint64_t id_absolute = current_batch_starting_sequence_id_;
  SingleSequence *sequence = NULL;
  std::string line;

  std::string sam_header = "";

  while (!ReadGZLine_(gzip_fp_, line)) {
    if (line.size() == 0) continue;

//    printf ("%s\n", line.c_str());

    if (line[0] == '@') {
      sam_header += line;

    } else {
      sequence = new SingleSequence();
      std::istringstream ss(line);
      std::string header;
      std::string seq;
      std::string qual;
      std::string cigar_string;

      uint32_t flag;
      std::string rname;
      int64_t pos;
      int32_t mapq;
      std::string rnext;
      int64_t pnext;
      int64_t tlen;

//      ss >> header >> aln.flag >> aln.rname >> aln.pos >> aln.mapq >> cigar_string >> aln.rnext >> aln.pnext >> aln.tlen >> seq >> qual;
      ss >> header >> flag >> rname >> pos >> mapq >> cigar_string >> rnext >> pnext >> tlen >> seq >> qual;

//      SequenceAlignment::SplitCigar(cigar_string, aln.cigar);

      // Load optional parameters from a SAM line.
      std::vector<std::string> optional;
      std::string opt_par;
      while (ss >> opt_par) {
        optional.push_back(opt_par);
      }
//      aln.set_optional(optional);

      SequenceAlignment aln(flag, rname, pos, mapq, cigar_string, rnext, pnext, tlen, optional);
      aln.ProcessOptional();

      if (qual == "*") {
        sequence->InitHeaderAndDataFromAscii((char *) header.c_str(),
                                             header.length(),
                                                  (int8_t *) seq.c_str(), seq.length(), id, id_absolute);
      } else {
        sequence->InitAllFromAscii((char *) header.c_str(), header.length(),
                                   (int8_t *) seq.c_str(), (int8_t *) qual.c_str(), seq.length(),
                                   id, id_absolute);
      }

      sequence->InitAlignment(aln);
      if (randomize_non_acgt_bases == true) {
        sequence->RandomizeNonACGTBases();
      }
      if (convert_to_uppercase == true) {
        sequence->BasesToUppercase();
      }

      AddSequence(sequence, true);
      id += 1;  // Increment the relative sequence id counter.
      id_absolute += 1;

      if (num_seqs_to_load > 0 && id >= num_seqs_to_load)  // Batch loading stopping condition.
        break;
      if (megabytes_to_load > 0 && ConvertFromBytes(MEMORY_UNIT_MEGABYTE, current_data_size_) >= megabytes_to_load)  // Batch loading stopping condition.
        break;
    }
  }

  return id;
}

void SequenceFile::Sort() {
  std::sort(sequences_.begin(), sequences_.end(), seq_sort_key());
}
