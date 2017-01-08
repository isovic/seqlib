#include "edlib.h"
#include "edlibcigar.h"
#include <algorithm>
#include <string>
#include <vector>

int edlibAlignmentToCigar(unsigned char* alignment, int alignmentLength,
                          int cigarFormat, char** cigar_) {
    *cigar_ = NULL;
    if (cigarFormat != EDLIB_CIGAR_EXTENDED && cigarFormat != EDLIB_CIGAR_STANDARD) {
        return EDLIB_STATUS_ERROR;
    }

    // Maps move code from alignment to char in cigar.
    //                        0    1    2    3
    char moveCodeToChar[] = {'=', 'I', 'D', 'X'};
    if (cigarFormat == EDLIB_CIGAR_STANDARD) {
        moveCodeToChar[0] = moveCodeToChar[3] = 'M';
    }

    vector<char>* cigar = new vector<char>();
    char lastMove = 0;  // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
    for (int i = 0; i <= alignmentLength; i++) {
        // if new sequence of same moves started
        if (i == alignmentLength || (moveCodeToChar[alignment[i]] != lastMove && lastMove != 0)) {
            // Write number of moves to cigar string.
            int numDigits = 0;
            for (; numOfSameMoves; numOfSameMoves /= 10) {
                cigar->push_back('0' + numOfSameMoves % 10);
                numDigits++;
            }
            reverse(cigar->end() - numDigits, cigar->end());
            // Write code of move to cigar string.
            cigar->push_back(lastMove);
            // If not at the end, start new sequence of moves.
            if (i < alignmentLength) {
                // Check if alignment has valid values.
                if (alignment[i] > 3) {
                    delete cigar;
                    return EDLIB_STATUS_ERROR;
                }
                numOfSameMoves = 0;
            }
        }
        if (i < alignmentLength) {
            lastMove = moveCodeToChar[alignment[i]];
            numOfSameMoves++;
        }
    }
    cigar->push_back(0);  // Null character termination.
    *cigar_ = (char*) malloc(cigar->size() * sizeof(char));
    memcpy(*cigar_, &(*cigar)[0], cigar->size() * sizeof(char));
    delete cigar;

    return EDLIB_STATUS_OK;
}

int edlibAlignmentToCigar(unsigned char* alignment, int alignmentLength,
                          int cigarFormat, std::string &ret_cigar) {
    if (cigarFormat != EDLIB_CIGAR_EXTENDED && cigarFormat != EDLIB_CIGAR_STANDARD) {
        return EDLIB_STATUS_ERROR;
    }

    // Maps move code from alignment to char in cigar.
    //                        0    1    2    3
    char moveCodeToChar[] = {'=', 'I', 'D', 'X'};
    if (cigarFormat == EDLIB_CIGAR_STANDARD) {
        moveCodeToChar[0] = moveCodeToChar[3] = 'M';
    }

//    vector<char>* cigar = new vector<char>();
    ret_cigar = std::string("");
    char lastMove = 0;  // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
    for (int i = 0; i <= alignmentLength; i++) {
        // if new sequence of same moves started
        if (i == alignmentLength || (moveCodeToChar[alignment[i]] != lastMove && lastMove != 0)) {
            // Write number of moves to cigar string.
            int numDigits = 0;
            for (; numOfSameMoves; numOfSameMoves /= 10) {
                ret_cigar.push_back('0' + numOfSameMoves % 10);
                numDigits++;
            }
            std::reverse(ret_cigar.end() - numDigits, ret_cigar.end());
            // Write code of move to cigar string.
            ret_cigar.push_back(lastMove);
            // If not at the end, start new sequence of moves.
            if (i < alignmentLength) {
                // Check if alignment has valid values.
                if (alignment[i] > 3) {
                    return EDLIB_STATUS_ERROR;
                }
                numOfSameMoves = 0;
            }
        }
        if (i < alignmentLength) {
            lastMove = moveCodeToChar[alignment[i]];
            numOfSameMoves++;
        }
    }
    ret_cigar.push_back(0);

    return EDLIB_STATUS_OK;
}

int edlibCalcEditDistance(
        const unsigned char* query, int queryLength,
        const unsigned char* target, int targetLength,
        int alphabetLength, int k, EdlibAlignMode edlib_mode,
        bool findStartLocations, bool findAlignment,
        int* bestScore, int** endLocations, int** startLocations, int* numLocations,
        unsigned char** alignment, int* alignmentLength, int *ret_k) {

  EdlibAlignTask task = (findAlignment == true) ? EDLIB_TASK_PATH :
                        (findStartLocations == true && findAlignment == false) ? EDLIB_TASK_LOC :
                        EDLIB_TASK_DISTANCE;

  EdlibAlignResult result = edlibAlign((const char *) query, queryLength, (const char *) target, targetLength, edlibNewAlignConfig(k, edlib_mode, task));

  if (result.numLocations == 0) {
    edlibFreeAlignResult(result);
    edlibFreeAlignResult(result);
    return EDLIB_STATUS_ERROR;
  }

  *bestScore = result.editDistance;
  *numLocations = result.numLocations;

  *endLocations = NULL;
  *startLocations = NULL;
  *alignment = NULL;

  int *starts = NULL;
  if (findStartLocations == true) {
    starts = (int *) malloc(sizeof(int) * result.numLocations);
  }

  int *ends = (int *) malloc(sizeof(int) * result.numLocations);

  for (int i=0; i<result.numLocations; i++) {
    if (findStartLocations == true) {
      starts[i] = result.startLocations[i];
    }

    ends[i] = result.endLocations[i];
  }

  *alignmentLength = result.alignmentLength;
  *startLocations = starts;
  *endLocations = ends;

  if (findAlignment == true) {
    *alignment = (unsigned char *) malloc(sizeof(unsigned char) * (result.alignmentLength));
    memcpy(*alignment, result.alignment, sizeof(char) * result.alignmentLength);
  }

  if (ret_k != NULL) {
    *ret_k = -1;
  }

  edlibFreeAlignResult(result);

  return 0;
}

