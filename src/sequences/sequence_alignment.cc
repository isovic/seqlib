/*
 * sequence_alignment.cc
 *
 *  Created on: Feb 14, 2016
 *      Author: isovic
 */

#include "sequence_alignment.h"
#include <math.h>

SequenceAlignment::SequenceAlignment() {
  qname = rname = cigar = rnext = "";
  flag = 4;
  pos = pnext = tlen = as = 0;
  mapq = 0;
  evalue = 0.0;
  optional.clear();
}

SequenceAlignment::~SequenceAlignment() {

}

void SequenceAlignment::CopyFrom(const SequenceAlignment& aln) {
  this->qname = aln.qname;
  this->flag = aln.flag;
  this->rname = aln.rname;
  this->pos = aln.pos;
  this->mapq = aln.mapq;
  this->cigar = aln.cigar;
  this->rnext = aln.rnext;
  this->pnext = aln.pnext;
  this->tlen = aln.tlen;
  this->as = aln.as;
  this->evalue = aln.evalue;
  this->optional = aln.optional;
}

void SequenceAlignment::ProcessOptional() {
  for (int32_t i=0; i<optional.size(); i++) {
    // Optional parameters are split with ':', but on specific characters.
//    if (optional[i].substr(0, 2) == "AS") {
//    }
  }
}

int SequenceAlignment::GetSplitCigar(std::vector<CigarOp>& ret) const {
  SplitCigar(cigar, ret);
  return 0;
}

int SequenceAlignment::SplitCigar(const std::string &cigar_str, std::vector<CigarOp>& ret) {
  ret.clear();
  CigarOp op;
  int32_t digit_count = 0;
  const char *first_num = NULL;
  for (int64_t i=0; i<cigar_str.size(); i++) {
    if (isalpha(cigar_str[i]) || cigar_str[i] == '=') {
      op.op = cigar_str[i];
      sscanf(first_num, "%d", &op.count);
      ret.push_back(op);
      first_num = NULL;
    } else if (first_num == NULL) {
      first_num = &(cigar_str[i]);
    }
  }
  return 0;
}

int SequenceAlignment::CalcReferenceLengthFromCigar(const std::vector<CigarOp>& split_cigar, int64_t& ret_ref_len) {
  int64_t len = 0;
  for (int64_t i=0; i<split_cigar.size(); i++) {
    if (is_cigar_ref(split_cigar[i].op)) { len += split_cigar[i].count; }
  }
  ret_ref_len = len;
  return 0;
}

bool SequenceAlignment::IsMapped() const {
  return (!(flag & 4));
}

int SequenceAlignment::CalcQueryLengthFromCigar(const std::vector<CigarOp>& split_cigar, int64_t& ret_query_len) {
  int64_t len = 0;
  for (int64_t i=0; i<split_cigar.size(); i++) {
    if (is_cigar_read(split_cigar[i].op)) { len += split_cigar[i].count; }
  }
  ret_query_len = len;
  return 0;
}


int64_t SequenceAlignment::GetReferenceLengthFromCigar() const {
  int64_t len = 0;
  std::vector<CigarOp> split_cigar;
  SplitCigar(cigar, split_cigar);
  CalcReferenceLengthFromCigar(split_cigar, len);
  return len;
}

int64_t SequenceAlignment::GetQueryLengthFromCigar() const {
  int64_t len = 0;
  std::vector<CigarOp> split_cigar;
  SplitCigar(cigar, split_cigar);
  CalcQueryLengthFromCigar(split_cigar, len);
  return len;
}
