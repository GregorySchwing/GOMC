#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "MultiParticle.h"


class MPTest : public ::testing::Test {
 protected:
  void SetUp() override {
     lb_new = XYZ(1.0, 1.0, 1.0);
     lb_old = XYZ(1.0, 1.0, 1.0);
     k = XYZ(1.0, 1.0, 1.0);
     max = 0.0;
  }

  // void TearDown() override {}
    XYZ lb_new; 
    XYZ lb_old;
    XYZ k; 
    double max;
    //StaticVals statV;
    //System sys;
    MultiParticle * mv;

};

TEST_F(MPTest, IsInitializedToOne) {
  EXPECT_EQ(lb_new.x, 1.0);
}

TEST_F(MPTest, CalcWRatio) {
  EXPECT_EQ(mv->CalculateWRatio(lb_new, lb_old, k, max), 1.0);
}

// Simple test, does not use gmock
/*TEST(Dummy, foobar)
{
    EXPECT_EQ(MultiParticle::CalculateWRatio(), 0);
}
*/

// Real class we want to mock
class TeaBreak
{
public:
    virtual ~TeaBreak() {}

    // Return minutes taken to make the drinks
    int morningTea()
    {
        return makeCoffee(true,  1) +
               makeCoffee(false, 0.5) +
               makeHerbalTea();
    }

private:
    virtual int makeCoffee(bool milk, double sugars) = 0;
    virtual int makeHerbalTea() = 0;
};

// Mock class
class MockTeaBreak : public TeaBreak
{
public:
    MOCK_METHOD2(makeCoffee,    int(bool milk, double sugars));
    MOCK_METHOD0(makeHerbalTea, int());
};


using ::testing::Return;
using ::testing::_;

// Mocked test
TEST(TeaBreakTest, MorningTea)
{
    MockTeaBreak  teaBreak;
    EXPECT_CALL(teaBreak, makeCoffee(_,_))
        .WillOnce(Return(2))
        .WillOnce(Return(1));
    EXPECT_CALL(teaBreak, makeHerbalTea())
        .WillOnce(Return(3));

    EXPECT_LE(teaBreak.morningTea(), 6);
}
