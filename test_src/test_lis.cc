#include <stdint.h>
#include <vector>
#include <assert.h>

#include "gtest/gtest.h"

#include "../src/lis.hpp"

TEST(LIS, LisTest1) {
  std::vector<int> v{};
  auto lis = LIS<int32_t>(v);
  std::vector<int32_t> expected = {};
  EXPECT_EQ(lis, expected);
}

TEST(LIS, LisTest2) {
  std::vector<int> v{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  auto lis = LIS<int32_t>(v);
  std::vector<int32_t> expected = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  EXPECT_EQ(lis, expected);
}

TEST(LIS, LisTest3) {
  std::vector<int> v{2, 5, 3, 7, 11, 8, 10, 13, 6};
  auto lis = LIS<int32_t>(v);
  std::vector<int32_t> expected = {2, 3, 7, 8, 10, 13};
  EXPECT_EQ(lis, expected);
}

TEST(LIS, LisTest4) {
  std::vector<int> v{1};
  auto lis = LIS<int32_t>(v);
  std::vector<int32_t> expected = {1};
  EXPECT_EQ(lis, expected);
}

TEST(LIS, LisTest5) {
  std::vector<int> v{3,1,4,1,5,9,2,6,5,3,5,8,9,7,9};
  auto lis = LIS<int32_t>(v);
  std::vector<int32_t> expected = {1,2,3,5,8,9};
  EXPECT_EQ(lis, expected);
}

class Hit {
 public:
  Hit(int32_t _x, int32_t _y) : x(_x), y(_y) { }
  int32_t x, y;
  bool operator==(const Hit& op) const {  // Needed for EXPECT_EQ.
    return (this->x == op.x && this->y == op.y);
  }
};

TEST(LIS, LisTest6) {
  // Test a custom data structure and comparison operator.
  std::vector<Hit> v = {Hit(1, 3),
                        Hit(2, 1),
                        Hit(3, 4),
                        Hit(4, 1),
                        Hit(5, 5),
                        Hit(6, 9),
                        Hit(7, 2),
                        Hit(8, 6),
                        Hit(9, 5),
                        Hit(10, 3),
                        Hit(11, 5),
                        Hit(12, 8),
                        Hit(13, 9),
                        Hit(14, 7),
                        Hit(15, 9)
                      };

  // Define a custom comparison function.
  auto lis = LIS<Hit>(v, [](const Hit& a, const Hit& b) { return a.y < b.y; });

  std::vector<Hit> expected = {
                                Hit(2, 1),
                                Hit(7, 2),
                                Hit(10, 3),
                                Hit(11, 5),
                                Hit(12, 8),
                                Hit(13, 9)
                              };

  EXPECT_EQ(lis, expected);
}
