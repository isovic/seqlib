/*
 * sequence_alignment.cc
 *
 *  Created on: Feb 14, 2016
 *      Author: isovic
 */

#include "sequence_alignment.h"
#include <math.h>
#include <sstream>

SequenceAlignment::SequenceAlignment() {
//  qname = "";
  rname = rnext = "";
  flag = 4;
  pos = pnext = tlen = as = 0;
  mapq = 0;
  evalue = 0.0;
  optional.clear();
  cigar.clear();
}

SequenceAlignment::~SequenceAlignment() {

}

void SequenceAlignment::CopyFrom(const SequenceAlignment& aln) {
//  this->qname = aln.qname;
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
//
//int SequenceAlignment::GetSplitCigar(std::vector<CigarOp>& ret) const {
//  SplitCigar(cigar, ret);
//  return 0;
//}

int SequenceAlignment::SplitCigar(const std::string &cigar_str, std::vector<CigarOp>& ret) {
  ret.clear();
  CigarOp op;
  int32_t digit_count = 0;
  int64_t pos_ref = 0, pos_query = 0;
  const char *first_digit = NULL;
  for (int64_t i=0; i<cigar_str.size(); i++) {
    if (isalpha(cigar_str[i]) || cigar_str[i] == '=') {
      op.pos_ref = pos_ref;
      op.pos_query = pos_query;
      op.op = cigar_str[i];
      sscanf(first_digit, "%d", &op.count);
      ret.push_back(op);
      first_digit = NULL;
      if (is_cigar_ref(op.op)) pos_ref += op.count;
      if (is_cigar_read(op.op)) pos_query += op.count;
    } else if (first_digit == NULL) {
      first_digit = &(cigar_str[i]);
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

int64_t SequenceAlignment::FindBasePositionOnRead(const std::vector<CigarOp>&  split_cigar, int64_t pos_on_ref, int64_t *cigar_id) const {
  int64_t aligned_pos = (this->pos > 0) ? (this->pos - 1) : 0;      // pos field of the SAM file is 1-based. If the value is <= 0, then something is not set right, and ignore this value.
  if (split_cigar.size() > 0 && (pos_on_ref < (split_cigar[0].pos_ref + aligned_pos))) {
    return -1;
  }

  int64_t i=0;
  for (i=0; i<split_cigar.size(); i++) {
    char op = split_cigar[i].op;
    int64_t cig_pos_ref = split_cigar[i].pos_ref + aligned_pos;
    if (is_cigar_match(op) || is_cigar_del(op)) {
      if (pos_on_ref >= cig_pos_ref && pos_on_ref < (cig_pos_ref + split_cigar[i].count)) {
        if (cigar_id != NULL) *cigar_id = i;
        if (is_cigar_match(op)) { return (split_cigar[i].pos_query + (pos_on_ref - cig_pos_ref)); }
        else { return (split_cigar[i].pos_query); }
      }
    }
  }
//  int64_t ref_len;
//  CalcReferenceLengthFromCigar(split_cigar, ref_len);
//  fprintf (stderr, "WARNING: i = %ld, split_cigar.size() = %ld, pos_on_ref = %ld, pos = %ld, length_on_ref = %ld, end_on_ref = %ld\n", i, split_cigar.size(), pos_on_ref, (pos - 1), ref_len, (pos - 1 + ref_len - 1));
  return -2;
}

int64_t SequenceAlignment::FindBasePositionOnRef(const std::vector<CigarOp>& split_cigar, int64_t pos_on_read, int64_t *cigar_id) const {
  int64_t aligned_pos = (this->pos > 0) ? (this->pos - 1) : 0;      // pos field of the SAM file is 1-based. If the value is <= 0, then something is not set right, and ignore this value.
  if (pos_on_read < 0) {
    return -1;
  }

  for (int64_t i=0; i<split_cigar.size(); i++) {
    char op = split_cigar[i].op;
//    int64_t cig_pos_ref = split_cigar[i].pos_ref + aligned_pos;
    if (is_cigar_match(op)) {
      if (pos_on_read >= split_cigar[i].pos_query && pos_on_read < (split_cigar[i].pos_query + split_cigar[i].count)) {
        if (cigar_id != NULL) *cigar_id = i;
        return (split_cigar[i].pos_ref + aligned_pos + (pos_on_read - split_cigar[i].pos_query));
      }
    } else if (is_cigar_ins(op)) {
      if (pos_on_read >= split_cigar[i].pos_query && pos_on_read < (split_cigar[i].pos_query + split_cigar[i].count)) {
        if (cigar_id != NULL) *cigar_id = i;
        return (split_cigar[i].pos_ref + aligned_pos);
      }
    }
  }
  return -2;
}

std::string SequenceAlignment::GetCigarString() const {
  return MakeCigarString(cigar);
}

std::string SequenceAlignment::MakeCigarString(const std::vector<CigarOp>& split_cigar) {
  std::stringstream ss;
  for (int64_t i=0; i<split_cigar.size(); i++) {
    ss << split_cigar[i].count << split_cigar[i].op;
  }
  return ss.str();
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
  CalcReferenceLengthFromCigar(cigar, len);
  return len;
}

int64_t SequenceAlignment::GetQueryLengthFromCigar() const {
  int64_t len = 0;
  CalcQueryLengthFromCigar(cigar, len);
  return len;
}

int64_t SequenceAlignment::GetClippedBasesFront() const {
  int64_t len = 0;
  CalcClippedBasesFront(cigar, len);
  return len;
}

int64_t SequenceAlignment::GetClippedBasesBack() const {
  int64_t len = 0;
  CalcClippedBasesBack(cigar, len);
  return len;
}

int SequenceAlignment::CalcClippedBasesFront(const std::vector<CigarOp>& split_cigar, int64_t& ret_clip_len) {
  int64_t clip_len = 0;
  for (int64_t i=0; i<split_cigar.size(); i++) {
    if (is_cigar_soft(split_cigar[i].op)) { clip_len += split_cigar[i].count; }
    else if (is_cigar_hard(split_cigar[i].op)) { }
    else { break; }
  }
  ret_clip_len = clip_len;
  return 0;
}

int SequenceAlignment::CalcClippedBasesBack(const std::vector<CigarOp>& split_cigar, int64_t& ret_clip_len) {
  int64_t clip_len = 0;
  for (int64_t i=(split_cigar.size()-1); i>=0; i--) {
    if (is_cigar_soft(split_cigar[i].op)) { clip_len += split_cigar[i].count; }
    else if (is_cigar_hard(split_cigar[i].op)) { }
    else { break; }
  }
  ret_clip_len = clip_len;
  return 0;
}
