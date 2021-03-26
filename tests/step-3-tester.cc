#include <gtest/gtest.h>

#include <fstream>

#include "step-3.h"

using namespace dealii;

// Test Fixture for step-3
class Step3Tester : public ::testing::Test, public Step3
{
public:
  Step3Tester() = default;
};


TEST_F(Step3Tester, Ex3)
{
  make_grid();
}

TEST_F(Step3Tester, Ex4)
{
  modify_bdary_cond = true;
  run();
}

TEST_F(Step3Tester, Ex5)
{
  modify_bdary_cond = true;
  modify_bdary_data = true;
  run();
}

TEST_F(Step3Tester, Ex6)
{
  source_term = 0.0; // set source term f = 0
  run();
}

TEST_F(Step3Tester, Ex7)
{
  l_shaped = true; // solve over L-shaped domain
  run();
}