#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "MultiParticle.h"
#include <limits>

class MPTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Basic Example
    lb_new = XYZ(1.0, 1.0, 1.0);
    lb_old = XYZ(1.0, 1.0, 1.0);
    k = XYZ(1.0, 1.0, 1.0);
    max = 0.0;

    // Tiny Max, Giant Old,
    double giantOLD = 1E296;
    lb_new_TMGO = XYZ(   std::numeric_limits<double>::infinity(), 
                    std::numeric_limits<double>::infinity(), 
                    std::numeric_limits<double>::infinity());
    giant_lb_old = XYZ(giantOLD, giantOLD, giantOLD);
    k = XYZ(1.0, 1.0, 1.0);
    tinyMax = 1E-307;


  // Giant Max, Tiny Old,  INF New
    double tinyOLD = 1E-307;
    lb_new_GMTO = XYZ(   std::numeric_limits<double>::infinity(), 
                    std::numeric_limits<double>::infinity(), 
                    std::numeric_limits<double>::infinity());
    tiny_lb_old = XYZ(tinyOLD, tinyOLD, tinyOLD);
    k = XYZ(1.0, 1.0, 1.0);
    giantMax = 1E296;
  }

  // void TearDown() override {}
  // Basic Example
    XYZ lb_new; 
    XYZ lb_old;
    XYZ k; 
    double max;

  // Tiny Max, Giant Old
    XYZ lb_new_TMGO; 
    XYZ giant_lb_old;
    XYZ k_TMGO; 
    double tinyMax;


  // Giant Max, Tiny Old
    XYZ lb_new_GMTO; 
    XYZ tiny_lb_old;
    XYZ lb_new_inf; 
    XYZ k_GMTO; 
    double giantMax;

    MultiParticle * mv;

};

TEST_F(MPTest, IsInitializedCorrectly) {
  EXPECT_EQ(lb_new_TMGO.x, std::numeric_limits<double>::infinity());
  EXPECT_EQ(lb_new_TMGO.x * 0, 0);

  EXPECT_GT(giant_lb_old.x * tinyMax, MIN_FORCE);
  EXPECT_LT(giant_lb_old.x * tinyMax, MAX_FORCE);
  EXPECT_GT(tiny_lb_old.x * giantMax, MIN_FORCE);
  EXPECT_LT(tiny_lb_old.x * giantMax, MAX_FORCE);
}

TEST_F(MPTest, CalcWRatio) {
  EXPECT_EQ(mv->CalculateWRatio(lb_new, lb_old, k, max), 1.0);
    // Tiny Max, Giant Old
  EXPECT_NE(mv->CalculateWRatio(lb_new_TMGO, giant_lb_old, k, tinyMax), std::numeric_limits<double>::infinity());
    // Giant Max, Tiny Old
  EXPECT_NE(mv->CalculateWRatio(lb_new_GMTO, tiny_lb_old, k, giantMax), std::numeric_limits<double>::infinity());

  ASSERT_TRUE(isfinite(mv->CalculateWRatio(lb_new_GMTO, tiny_lb_old, k, giantMax)));

}