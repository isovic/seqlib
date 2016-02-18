/*
 * sequence_alignment_test.cc
 *
 *  Created on: Feb 15, 2016
 *      Author: isovic
 */

#include "sequence_alignment.h"
#include <assert.h>

void TEST_SEQUENCE_ALIGNMENT() {
  SequenceAlignment aln;

  aln.cigar = "12M1I10=2X10=25S";

  std::vector<CigarOp> split_ops;
  aln.GetSplitCigar(split_ops);

  std::vector<CigarOp> test_ops;
  CigarOp op;
  op.op = 'M';  op.count = 12;  test_ops.push_back(op);
  op.op = 'I';  op.count = 1;   test_ops.push_back(op);
  op.op = '=';  op.count = 10;  test_ops.push_back(op);
  op.op = 'X';  op.count = 2;  test_ops.push_back(op);
  op.op = '=';  op.count = 10;  test_ops.push_back(op);
  op.op = 'S';  op.count = 25;  test_ops.push_back(op);

  assert(split_ops.size() == test_ops.size());
  for (int32_t i=0; i<test_ops.size(); i++) {
    assert(split_ops[i].op == test_ops[i].op);
    assert(split_ops[i].count == test_ops[i].count);
  }

  fprintf (stderr, "%s passed.\n", __FUNCTION__);
}
