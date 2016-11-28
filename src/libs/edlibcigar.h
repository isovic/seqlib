#ifndef EDLIBCIGAR_HEADER_H
#define EDLIBCIGAR_HEADER_H

#include <stdint.h>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <cstring>
#include <cassert>
#include <stdio.h>

using namespace std;

typedef uint64_t Word;
static const int WORD_SIZE = sizeof(Word) * 8; // Size of Word in bits
static const Word WORD_1 = (Word)1;
static const Word HIGH_BIT_MASK = WORD_1 << (WORD_SIZE - 1);  // 100..00

/**
 * Builds cigar string from given alignment sequence.
 * @param [in] alignment  Alignment sequence.
 *     0 stands for match.
 *     1 stands for insertion to target.
 *     2 stands for insertion to query.
 *     3 stands for mismatch.
 * @param [in] alignmentLength
 * @param [in] cigarFormat
 *     If EDLIB_CIGAR_EXTENDED, extended cigar is returned.
 *     If EDLIB_CIGAR_STANDARD, standard cigar is returned (contains only I, D and M).
 * @param [out] cigar  Will contain cigar string.
 *     I stands for insertion.
 *     D stands for deletion.
 *     X stands for mismatch. (used only in extended format)
 *     = stands for match. (used only in extended format)
 *     M stands for (mis)match. (used only in standard format)
 *     String is null terminated.
 *     Needed memory is allocated and given pointer is set to it.
 *     Do not forget to free it later using free()!
 * @return Status code.
 */
int edlibAlignmentToCigar(unsigned char* alignment, int alignmentLength,
                          int cigarFormat, char** cigar);

int edlibCalcEditDistance(
        const unsigned char* query, int queryLength,
        const unsigned char* target, int targetLength,
        int alphabetLength, int k, EdlibAlignMode edlib_mode,
        bool findStartLocations, bool findAlignment,
        int* bestScore, int** endLocations, int** startLocations, int* numLocations,
        unsigned char** alignment, int* alignmentLength, int *ret_k);

#endif
