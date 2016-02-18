/*
 * sequence_alignment.h
 *
 *  Created on: Feb 14, 2016
 *      Author: isovic
 */

#ifndef CODEBASE_SEQLIB_SRC_SEQUENCES_SEQUENCE_ALIGNMENT_H_
#define CODEBASE_SEQLIB_SRC_SEQUENCES_SEQUENCE_ALIGNMENT_H_

#include <stdint.h>
#include <string>
#include <vector>

#define is_cigar_op(x)  (x == 'M' || x == '=' || x == 'X' || x == 'I' || x == 'D' || x == 'S' || x == 'H')
#define is_cigar_match(x)  (x == 'M' || x == '=' || x == 'X')
#define is_cigar_ins(x)  (x == 'I')
#define is_cigar_del(x)  (x == 'D')
#define is_cigar_soft(x)  (x == 'S')
#define is_cigar_hard(x)  (x == 'H')
#define is_cigar_ref(x)  (x == 'M' || x == '=' || x == 'X' || x == 'D')
#define is_cigar_read(x)  (x == 'M' || x == '=' || x == 'X' || x == 'I' || x == 'S')

typedef struct {
  char op = '-';
  int32_t count = 0;
} CigarOp ;

class SequenceAlignment {
 public:
  SequenceAlignment();
  ~SequenceAlignment();

  void CopyFrom(const SequenceAlignment &aln);
  // Parses optional parameters (from member optional) to obtain values for 'as', 'evalue', etc.
  void ProcessOptional();
  int GetSplitCigar(std::vector<CigarOp> &ret) const;
  int64_t GetReferenceLengthFromCigar() const;
  int64_t GetQueryLengthFromCigar() const;
  bool IsMapped() const;

  static int SplitCigar(const std::string &cigar_str, std::vector<CigarOp>& ret);

  static int CalcReferenceLengthFromCigar(const std::vector<CigarOp>& split_cigar, int64_t &ret_ref_len);
  static int CalcQueryLengthFromCigar(const std::vector<CigarOp>& split_cigar, int64_t &ret_query_len);

  std::string qname;    // Field #1.
  uint32_t flag;        // Field #2.
  std::string rname;    // Field #3.
  int64_t pos;          // Field #4.
  int32_t mapq;         // Field #5.
  std::string cigar;    // Field #6.
  std::string rnext;    // Field #7.
  int64_t pnext;        // Field #8.
  int64_t tlen;         // Field #9.
  // std::string seq;   // Field #10. Skipped because SingleSequence holds this info.
  // std::string qual;  // Field #11. Skipped because SingleSequence holds this info.

  // Optional fields in the SAM format:
  int64_t as;           // Alignment score.
  double evalue;        // E-value. There is no dedicated field in the SAM format, but GraphMap uses ZE to specify the E-value.
  std::vector<std::string> optional;  // Raw values (strings) of optional fields, not explicitly converted to expected values;
};

#endif /* CODEBASE_SEQLIB_SRC_SEQUENCES_SEQUENCE_ALIGNMENT_H_ */
