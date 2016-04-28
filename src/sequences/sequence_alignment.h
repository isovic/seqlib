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
  int64_t pos_ref = -1;
  int64_t pos_query = - 1;
} CigarOp ;

class SequenceAlignment {
 public:
  SequenceAlignment();
  SequenceAlignment(uint32_t flag, std::string &rname, int64_t pos, int32_t mapq, std::string &cigar_string, std::string &rnext, int64_t pnext, int64_t tlen, std::vector<std::string> &optional);
  ~SequenceAlignment();

  void CopyFrom(const SequenceAlignment &aln);
  // Parses optional parameters (from member optional) to obtain values for 'as', 'evalue', etc.
  void ProcessOptional();
  int64_t GetReferenceLengthFromCigar() const;
  int64_t GetQueryLengthFromCigar() const;
  // Returns the position if found. If the requested position begins before the position of the first base, the function returns -1. If the requested position is behind the last base, the function returns -2.
  // Parameter cigar_id is the ID (ordinal number) of the cigar operation where the position was found in. Optional.
  int64_t FindBasePositionOnRead(int64_t pos, int64_t *cigar_id=NULL) const;
  // Returns the position if found. If the requested position begins before the position of the first base, the function returns -1. If the requested position is behind the last base, the function returns -2.
  // Parameter cigar_id is the ID (ordinal number) of the cigar operation where the position was found in. Optional.
  int64_t FindBasePositionOnRef(int64_t pos, int64_t *cigar_id=NULL) const;
  std::string GetCigarString() const;
  void SetCigarFromString(std::string &cigar_str);

  bool IsMapped() const;

  static int CalcReferenceLengthFromCigar(const std::vector<CigarOp>& split_cigar, int64_t &ret_ref_len);
  static int CalcQueryLengthFromCigar(const std::vector<CigarOp>& split_cigar, int64_t &ret_query_len);
  static std::string MakeCigarString(const std::vector<CigarOp>& split_cigar);
  static int SplitCigar(const std::string &cigar_str, std::vector<CigarOp>& ret);

  int64_t get_as() const;
  void set_as(int64_t as);
  const std::vector<CigarOp>& get_cigar() const;
  void set_cigar(const std::vector<CigarOp>& cigar);
  double get_evalue() const;
  void set_evalue(double evalue);
  uint32_t get_flag() const;
  void set_flag(uint32_t flag);
  int32_t get_mapq() const;
  void set_mapq(int32_t mapq);
  const std::vector<std::string>& get_optional() const;
  void set_optional(const std::vector<std::string>& optional);
  int64_t get_pnext() const;
  void set_pnext(int64_t pnext);
  int64_t get_pos() const;
  void set_pos(int64_t pos);
  const std::string& get_rname() const;
  void set_rname(const std::string& rname);
  const std::string& get_rnext() const;
  void set_rnext(const std::string& rnext);
  int64_t get_tlen() const;
  void set_tlen(int64_t tlen);

 private:

//  std::string qname;    // Field #1.
  uint32_t flag_;        // Field #2.
  std::string rname_;    // Field #3.
  int64_t pos_;          // Field #4. At the moment, 1-based, as in SAM file.
  int32_t mapq_;         // Field #5.
//  std::string cigar;    // Field #6.
  std::vector<CigarOp> cigar_;
  std::string rnext_;    // Field #7.
  int64_t pnext_;        // Field #8.
  int64_t tlen_;         // Field #9.
  // std::string seq;   // Field #10. Skipped because SingleSequence holds this info.
  // std::string qual;  // Field #11. Skipped because SingleSequence holds this info.

  // Optional fields in the SAM format:
  int64_t as_;           // Alignment score.
  double evalue_;        // E-value. There is no dedicated field in the SAM format, but GraphMap uses ZE to specify the E-value.
  std::vector<std::string> optional_;  // Raw values (strings) of optional fields, not explicitly converted to expected values;
};

#endif /* CODEBASE_SEQLIB_SRC_SEQUENCES_SEQUENCE_ALIGNMENT_H_ */
