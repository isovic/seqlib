/*
 * sequence_alignment_test.cc
 *
 *  Created on: Feb 15, 2016
 *      Author: isovic
 */

#include "sequence_alignment.h"
#include <assert.h>

void TEST_CLASS_SEQUENCE_ALIGNMENT() {
  SequenceAlignment aln;

//  aln.cigar = "12M3I10=2X10=7D5=25S";
//  std::vector<CigarOp> split_ops;
//  aln.GetSplitCigar(split_ops);

  std::string cigar_string = "12M3I10=2X10=7D5=25S";
  std::vector<CigarOp> split_ops;
//  SequenceAlignment::SplitCigar(cigar_string, aln.get_cigar());
  SequenceAlignment::SplitCigar(cigar_string, split_ops);
  aln.SetCigarFromString(cigar_string);

  std::vector<CigarOp> test_ops;
  CigarOp op;
  op.op = 'M';  op.count = 12;  test_ops.push_back(op);
  op.op = 'I';  op.count = 3;   test_ops.push_back(op);
  op.op = '=';  op.count = 10;  test_ops.push_back(op);
  op.op = 'X';  op.count = 2;  test_ops.push_back(op);
  op.op = '=';  op.count = 10;  test_ops.push_back(op);
  op.op = 'D';  op.count = 7;  test_ops.push_back(op);
  op.op = '=';  op.count = 5;  test_ops.push_back(op);
  op.op = 'S';  op.count = 25;  test_ops.push_back(op);

  assert(split_ops.size() == test_ops.size());
  for (size_t i=0; i<test_ops.size(); i++) {
    assert(split_ops[i].op == test_ops[i].op);
    assert(split_ops[i].count == test_ops[i].count);
  }

  for (size_t i=0; i<split_ops.size(); i++) {
    printf ("%d %c, read = %lld, ref = %lld\n", split_ops[i].count, split_ops[i].op, split_ops[i].pos_query, split_ops[i].pos_ref);
  }

  int64_t ref_pos = 0;
  int64_t read_pos = 0;
  read_pos = 14;
  ref_pos = aln.FindBasePositionOnRef(read_pos);
  printf ("read_pos = %lld -> ref_pos = %lld\n", read_pos, ref_pos);
  assert(ref_pos == 12);
  read_pos = 27;
  ref_pos = aln.FindBasePositionOnRef(read_pos);
  printf ("read_pos = %lld -> ref_pos = %lld\n", read_pos, ref_pos);
  assert(ref_pos == 24);
  read_pos = 37;
  ref_pos = aln.FindBasePositionOnRef(read_pos);
  printf ("read_pos = %lld -> ref_pos = %lld\n", read_pos, ref_pos);
  assert(ref_pos == 41);

  ref_pos = 14;
  read_pos = aln.FindBasePositionOnRead(ref_pos);
  printf ("ref_pos = %lld -> read_pos = %lld\n", ref_pos, read_pos);
  assert(read_pos == 17);
  ref_pos = 27;
  read_pos = aln.FindBasePositionOnRead(ref_pos);
  printf ("ref_pos = %lld -> read_pos = %lld\n", ref_pos, read_pos);
  assert(read_pos == 30);
  ref_pos = 37;
  read_pos = aln.FindBasePositionOnRead(ref_pos);
  printf ("ref_pos = %lld -> read_pos = %lld\n", ref_pos, read_pos);
  assert(read_pos == 37);

//  fprintf (stderr, "FindBasePositionOnRef(%ld) = %ld\n", 37, ref_pos);

  fprintf (stderr, "%s passed.\n", __FUNCTION__);
}
