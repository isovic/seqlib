/*
 * sequence_gfa.cc
 *
 *  Created on: Apr 28, 2016
 *      Author: isovic
 */

#include <sstream>
#include <algorithm>

#include "sequence_gfa.h"

SequenceGFA::SequenceGFA() {
  // TODO Auto-generated constructor stub

}

SequenceGFA::~SequenceGFA() {
  // TODO Auto-generated destructor stub
}

const std::vector<GFAGoldenPath>& SequenceGFA::get_gpath() const {
  return gpath_;
}

void SequenceGFA::set_gpath(const std::vector<GFAGoldenPath>& gpath) {
  gpath_ = gpath;
}

const std::vector<GFAOverlap>& SequenceGFA::get_overlaps() const {
  return overlaps_;
}

void SequenceGFA::AddGoldenPath(GFAGoldenPath& gp) {
  gpath_.push_back(gp);
}

void SequenceGFA::AddOverlap(GFAOverlap& ov) {
  overlaps_.push_back(ov);
}

void SequenceGFA::ParseGoldenPath(const std::string &line, GFAGoldenPath &gpath) {
  std::istringstream ss(line);
  std::string keyword;
  std::string rbuffer;
  ss >> keyword >> gpath.utg_name >> gpath.utg_start >> rbuffer >> gpath.orientation >> gpath.inc_len;

  auto pos1 = rbuffer.find(':');
  if (pos1 == std::string::npos) { return; }
  gpath.read_name = rbuffer.substr(0, pos1);

  auto pos2 = rbuffer.find('-');
  if (pos2 == std::string::npos) { return; }
  std::string temp_start = rbuffer.substr((pos1+1), pos2);
  std::string temp_end = rbuffer.substr((pos2+1));

  std::istringstream ss_start(temp_start);
  ss_start >> gpath.read_start;
  std::istringstream ss_end(temp_end);
  ss_end >> gpath.read_end;
}

void SequenceGFA::ParseOverlap(const std::string &line, GFAOverlap &ov) {
  std::istringstream ss(line);
  std::string keyword;
  ss >> keyword >> ov.seg1_name >> ov.seg1_orient >> ov.seg2_name >> ov.seg2_orient;
}

void SequenceGFA::CopyFrom(const SequenceGFA& gfa) {
  gpath_ = gfa.gpath_;
  overlaps_ = gfa.overlaps_;
}

void SequenceGFA::set_overlaps(const std::vector<GFAOverlap>& overlaps) {
  overlaps_ = overlaps;
}
