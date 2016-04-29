/*
 * sequence_gfa.h
 *
 *  Created on: Apr 28, 2016
 *      Author: isovic
 */

#ifndef SEQUENCE_GFA_H_
#define SEQUENCE_GFA_H_

#include <string>
#include <vector>

// This is a copy/paste from the Miniasm Man page:
// ┌─────┬─────────────┬──────────────────────────────────────────────────────┐
// │Line │   Comment   │                     Fixed fields                     │
// ├─────┼─────────────┼──────────────────────────────────────────────────────┤
// │ H   │ Header      │ N/A                                                  │
// │ S   │ Segment     │ segName segSeq                                       │
// │ L   │ Overlap     │ segName1 segOri1 segName2 segOri2 ovlpCIGAR          │
// │ a   │ Golden path │ utgName utgStart readName:rStart-rEnd readOri incLen │
// └─────┴─────────────┴──────────────────────────────────────────────────────┘

// Description of the GFAGoldenPath as described in the Miniasm Man page:
// An `a' line indicates that the unitig subsequence in [utgStart,utgStart+incLen)
// is taken from read readName in region [rStart-1,rStart-1+incLen).  It is not a
// standard GFA line. An `x' line gives  a  brief  summary  of  each
//       unitig, which can be inferred from `S' and `a' lines.
class GFAGoldenPath {
 public:
  std::string utg_name;
  int64_t utg_start;
  std::string read_name;
  int64_t read_start;
  int64_t read_end;
  char orientation;
  int64_t inc_len;
};

class GFAOverlap {
 public:
  std::string seg1_name;
  char seg1_orient;
  std::string seg2_name;
  char seg2_orient;
};

class SequenceGFA {
 public:
  SequenceGFA();
  ~SequenceGFA();

  void AddGoldenPath(GFAGoldenPath& gp);
  void AddOverlap(GFAOverlap& ov);
  void CopyFrom(const SequenceGFA& gfa);

  static void ParseGoldenPath(const std::string &line, GFAGoldenPath &gpath);
  static void ParseOverlap(const std::string &line, GFAOverlap &ov);

  const std::vector<GFAGoldenPath>& get_gpath() const;
  void set_gpath(const std::vector<GFAGoldenPath>& gpath);
  const std::vector<GFAOverlap>& get_overlaps() const;
  void set_overlaps(const std::vector<GFAOverlap>& overlaps);

 private:
  std::vector<GFAGoldenPath> gpath_;
  std::vector<GFAOverlap> overlaps_;
};

#endif /* CODEBASE_SEQLIB_SRC_SEQUENCES_SEQUENCE_GFA_H_ */
