/*
 * TicToc.h
 *
 *  Created on: Nov 29, 2016
 *      Author: isovic
 */

#ifndef SRC_CONSENSUS_TicToc_H_
#define SRC_CONSENSUS_TicToc_H_

#include <time.h>

class TicToc {
 public:
  TicToc();
  ~TicToc();

  void start();
  void stop();
  double get_secs();
  double get_msecs();

 private:
  struct timespec start_, stop_;
};

TicToc::TicToc() {
}

TicToc::~TicToc() {
}

void TicToc::start() {
  clock_gettime(CLOCK_MONOTONIC_RAW, &start_);
}

void TicToc::stop() {
  clock_gettime(CLOCK_MONOTONIC_RAW, &stop_);
}

double TicToc::get_secs() {
  double delta_us = (stop_.tv_sec - start_.tv_sec) * 1000000.0 + (stop_.tv_nsec - start_.tv_nsec) / 1000.0;
  return delta_us / 1000000.0;
}

double TicToc::get_msecs() {
  double delta_us = (stop_.tv_sec - start_.tv_sec) * 1000000.0 + (stop_.tv_nsec - start_.tv_nsec) / 1000.0;
  return delta_us / 1000.0;
}

#endif /* SRC_CONSENSUS_TicToc_H_ */
