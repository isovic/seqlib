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
  rname_ = rnext_ = "";
  flag_ = 4;
  pos_ = pnext_ = tlen_ = as_ = 0;
  mapq_ = 0;
  evalue_ = 0.0;
  optional_.clear();
  cigar_.clear();
}

SequenceAlignment::SequenceAlignment(uint32_t flag, std::string &rname, int64_t pos, int32_t mapq, std::string &cigar_string, std::string &rnext, int64_t pnext, int64_t tlen, std::vector<std::string> &optional)
: flag_(flag), rname_(rname), pos_(pos), mapq_(mapq), rnext_(rnext), pnext_(pnext), tlen_(tlen), optional_(optional) {
  SplitCigar(cigar_string, cigar_);
  ProcessOptional();
}

SequenceAlignment::~SequenceAlignment() {

}

void SequenceAlignment::CopyFrom(const SequenceAlignment& aln) {
//  this->qname = aln.qname;
  this->flag_ = aln.flag_;
  this->rname_ = aln.rname_;
  this->pos_ = aln.pos_;
  this->mapq_ = aln.mapq_;
  this->cigar_ = aln.cigar_;
  this->rnext_ = aln.rnext_;
  this->pnext_ = aln.pnext_;
  this->tlen_ = aln.tlen_;
  this->as_ = aln.as_;
  this->evalue_ = aln.evalue_;
  this->optional_ = aln.optional_;
}

void SequenceAlignment::ProcessOptional() {
  for (int32_t i=0; i<optional_.size(); i++) {
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
  return (!(flag_ & 4));
}

int64_t SequenceAlignment::FindBasePositionOnRead(int64_t pos_on_ref, int64_t *cigar_id) const {
  int64_t aligned_pos = (this->pos_ > 0) ? (this->pos_ - 1) : 0;      // pos field of the SAM file is 1-based. If the value is <= 0, then something is not set right, and ignore this value.
  if (cigar_.size() > 0 && (pos_on_ref < (cigar_[0].pos_ref + aligned_pos))) {
    return -1;
  }

//  int64_t i=0;
//  for (i=0; i<cigar_.size(); i++) {
//    char op = cigar_[i].op;
//    int64_t cig_pos_ref = cigar_[i].pos_ref + aligned_pos;
//    if (is_cigar_match(op) || is_cigar_del(op)) {
//      if (pos_on_ref >= cig_pos_ref && pos_on_ref < (cig_pos_ref + cigar_[i].count)) {
//        if (cigar_id != NULL) *cigar_id = i;
//        if (is_cigar_match(op)) { return (cigar_[i].pos_query + (pos_on_ref - cig_pos_ref)); }
//        else { return (cigar_[i].pos_query); }
//      }
//    }
//  }
//
////  int64_t ref_len;
////  CalcReferenceLengthFromCigar(split_cigar, ref_len);
////  fprintf (stderr, "WARNING: i = %ld, split_cigar.size() = %ld, pos_on_ref = %ld, pos = %ld, length_on_ref = %ld, end_on_ref = %ld\n", i, split_cigar.size(), pos_on_ref, (pos - 1), ref_len, (pos - 1 + ref_len - 1));
//  return -2;

  int64_t left = 0;
  int64_t right = cigar_.size() - 1;
  int64_t found_op_id = -1;
//  printf ("\n");
  while (left <= right) {
    int64_t mid = (left + right) / 2;
//    printf ("  mid = %ld, first = %ld, last = %ld, cigar_[mid].pos_ref = %ld, cigar_[left].pos_ref = %ld, cigar_[right].pos_ref = %ld\n", mid, left, right, cigar_[mid].pos_ref, cigar_[left].pos_ref, cigar_[right].pos_ref);
    if (pos_on_ref < (cigar_[mid].pos_ref + aligned_pos)) {
      right = mid - 1;
//      printf ("    -> right = mid - 1\n");
    } else if (pos_on_ref >= (cigar_[mid].pos_ref + aligned_pos) && pos_on_ref < (cigar_[mid].pos_ref + aligned_pos + cigar_[mid].count)) { // Note, this might find insertions as well! Insertions are not bases on the reference, we need to handle this case below.
//    } else if (pos_on_ref >= (cigar_[mid].pos_ref + aligned_pos) &&
//        ((cigar_[mid].op == 'I' && pos_on_ref == (cigar_[mid].pos_ref + aligned_pos)) ||
//          (is_cigar_ref(cigar_[mid].op) && pos_on_ref < (cigar_[mid].pos_ref + aligned_pos + cigar_[mid].count)))) {
//      printf ("found! pos_on_ref = %ld, cigar_[mid].pos_ref + aligned_pos = %ld, (cigar_[mid].pos_ref + aligned_pos + cigar_[mid].count) = %ld\n", pos_on_ref, cigar_[mid].pos_ref + aligned_pos, (cigar_[mid].pos_ref + aligned_pos + cigar_[mid].count));
      found_op_id = mid;
//      printf ("    -> found!\n");
      break;
    } else {
      left = mid + 1;
//      printf ("    -> left = mid + 1\n");
    }
  }
  if (left > right || found_op_id < 0) {
    // Not found.
//    printf ("\nTu sam 1!\nleft = %ld, right = %ld, found_op_id = %ld\npos_on_ref = %ld\n", left, right, found_op_id, pos_on_ref);
//    printf ("aligned_pos = %ld\n", aligned_pos);
//    for (int64_t i=0; i<cigar_.size(); i++) {
//      printf ("[%ld] %c %d, %ld, %ld\n", i, cigar_[i].op, cigar_[i].count, cigar_[i].pos_query, cigar_[i].pos_ref + aligned_pos);
//    }
//    fflush(stdout);
    return -2;
  }
  int64_t found_m_op = -1;
  int64_t ops_start = found_op_id, ops_end = found_op_id + 1;
  // Be conservative, and expand the cigar op range a bit to the left. This compensates for the insertions in the binary search.
  for (int64_t i=found_op_id; i>=0; i--) {
    int64_t cig_pos_ref = cigar_[i].pos_ref + aligned_pos;
//    if (pos_on_ref >= (cig_pos_ref) && pos_on_ref < (cig_pos_ref + cigar_[i].count)) { ops_start = i; }
//    else { break; }
    if (cig_pos_ref < pos_on_ref) { break; }
    ops_start = i;
  }
  // Be conservative, and expand the cigar op range a bit to the right. This compensates for the insertions in the binary search.
  for (int64_t i=(found_op_id+1); i<cigar_.size(); i++) {
    int64_t cig_pos_ref = cigar_[i].pos_ref + aligned_pos;
    if (cig_pos_ref > pos_on_ref) { break; }
    ops_end = (i + 1);
//    if (pos_on_ref >= (cig_pos_ref) && pos_on_ref < (cig_pos_ref + cigar_[i].count)) { ops_end = (i + 1); }
//    else { break; }

//    if ((cigar_[i].pos_ref + aligned_pos) == pos_on_ref) {
//      ops_end = (i + 1);
//      if (ops_start == 505) {
//      printf ("Tu sam 2!\n");
//      fflush(stdout);
//      printf ("ops_end = %ld\n", ops_end);
//      exit(1);
//      }
//    }
//    else { break; }
  }

  for (int64_t i=ops_start; i<ops_end; i++) {
    char op = cigar_[i].op;
    int64_t cig_pos_ref = cigar_[i].pos_ref + aligned_pos;
    if (is_cigar_match(op) || is_cigar_del(op)) {
      if (pos_on_ref >= cig_pos_ref && pos_on_ref < (cig_pos_ref + cigar_[i].count)) {
        if (cigar_id != NULL) *cigar_id = i;
        if (is_cigar_match(op)) { return (cigar_[i].pos_query + (pos_on_ref - cig_pos_ref)); }
        else { return (cigar_[i].pos_query); }
      }
    }
  }

//  printf ("Whaaat!?\n");
//  fflush(stdout);
//
//  printf ("ops_start = %ld, ops_end = %ld\n", ops_start, ops_end);
//
//  printf ("\nTu sam 1!\nleft = %ld, right = %ld, found_op_id = %ld\npos_on_ref = %ld\n", left, right, found_op_id, pos_on_ref);
//  printf ("aligned_pos = %ld\n", aligned_pos);
//  for (int64_t i=0; i<cigar_.size(); i++) {
//    printf ("[%ld] %c %d, %ld, %ld\n", i, cigar_[i].op, cigar_[i].count, cigar_[i].pos_query, cigar_[i].pos_ref + aligned_pos);
//  }
//  fflush(stdout);

  return -2;

//    if (cigar_[i].pos_ref == pos_on_ref &&
//        (is_cigar_match(cigar_[i].op)) &&
//        (found_m_op == -1 || i < found_m_op)) {
////      found_m_op = i;
//    }
//  }
}

int64_t SequenceAlignment::FindBasePositionOnRef(int64_t pos_on_read, int64_t *cigar_id) const {
  int64_t aligned_pos = (this->pos_ > 0) ? (this->pos_ - 1) : 0;      // pos field of the SAM file is 1-based. If the value is <= 0, then something is not set right, and ignore this value.
  if (pos_on_read < 0) {
    return -1;
  }

  for (int64_t i=0; i<cigar_.size(); i++) {
    char op = cigar_[i].op;
//    int64_t cig_pos_ref = split_cigar[i].pos_ref + aligned_pos;
    if (is_cigar_match(op)) {
      if (pos_on_read >= cigar_[i].pos_query && pos_on_read < (cigar_[i].pos_query + cigar_[i].count)) {
        if (cigar_id != NULL) *cigar_id = i;
        return (cigar_[i].pos_ref + aligned_pos + (pos_on_read - cigar_[i].pos_query));
      }
    } else if (is_cigar_ins(op)) {
      if (pos_on_read >= cigar_[i].pos_query && pos_on_read < (cigar_[i].pos_query + cigar_[i].count)) {
        if (cigar_id != NULL) *cigar_id = i;
        return (cigar_[i].pos_ref + aligned_pos);
      }
    }
  }
  return -2;
}

std::string SequenceAlignment::GetCigarString() const {
  return MakeCigarString(cigar_);
}

void SequenceAlignment::SetCigarFromString(std::string& cigar_str) {
  SplitCigar(cigar_str, cigar_);
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
  CalcReferenceLengthFromCigar(cigar_, len);
  return len;
}

int64_t SequenceAlignment::GetQueryLengthFromCigar() const {
  int64_t len = 0;
  CalcQueryLengthFromCigar(cigar_, len);
  return len;
}

int64_t SequenceAlignment::get_as() const {
  return as_;
}

void SequenceAlignment::set_as(int64_t as) {
  this->as_ = as;
}

const std::vector<CigarOp>& SequenceAlignment::get_cigar() const {
  return cigar_;
}

void SequenceAlignment::set_cigar(const std::vector<CigarOp>& cigar) {
  this->cigar_ = cigar;
}

std::vector<CigarOp>& SequenceAlignment::cigar() {
  return cigar_;
}

double SequenceAlignment::get_evalue() const {
  return evalue_;
}

void SequenceAlignment::set_evalue(double evalue) {
  this->evalue_ = evalue;
}

uint32_t SequenceAlignment::get_flag() const {
  return flag_;
}

void SequenceAlignment::set_flag(uint32_t flag) {
  this->flag_ = flag;
}

int32_t SequenceAlignment::get_mapq() const {
  return mapq_;
}

void SequenceAlignment::set_mapq(int32_t mapq) {
  this->mapq_ = mapq;
}

const std::vector<std::string>& SequenceAlignment::get_optional() const {
  return optional_;
}

void SequenceAlignment::set_optional(const std::vector<std::string>& optional) {
  this->optional_ = optional;
}

std::vector<std::string>& SequenceAlignment::optional() {
  return optional_;
}

int64_t SequenceAlignment::get_pnext() const {
  return pnext_;
}

void SequenceAlignment::set_pnext(int64_t pnext) {
  this->pnext_ = pnext;
}

int64_t SequenceAlignment::get_pos() const {
  return pos_;
}

void SequenceAlignment::set_pos(int64_t pos) {
  this->pos_ = pos;
}

const std::string& SequenceAlignment::get_rname() const {
  return rname_;
}

void SequenceAlignment::set_rname(const std::string& rname) {
  this->rname_ = rname;
}

const std::string& SequenceAlignment::get_rnext() const {
  return rnext_;
}

void SequenceAlignment::set_rnext(const std::string& rnext) {
  this->rnext_ = rnext;
}

int64_t SequenceAlignment::get_tlen() const {
  return tlen_;
}

void SequenceAlignment::set_tlen(int64_t tlen) {
  this->tlen_ = tlen;
}

int64_t SequenceAlignment::GetClippedBasesFront() const {
  int64_t len = 0;
  CalcClippedBasesFront(cigar_, len);
  return len;
}

int64_t SequenceAlignment::GetClippedBasesBack() const {
  int64_t len = 0;
  CalcClippedBasesBack(cigar_, len);
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
