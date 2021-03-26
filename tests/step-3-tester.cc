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


TEST_F(Step3Tester, MakeGrid)
{
  make_grid();
}

TEST_F(Step3Tester, ModifyBoundaryConditions)
{
  modify_bdary_cond = true;
  run();
}

TEST_F(Step3Tester, ModifyBoundaryData)
{
  modify_bdary_cond = true;
  modify_bdary_data = true;
  run();
}

TEST_F(Step3Tester, solve_laplace_equation)
{
  source_term = 0.0; // set source term f = 0
  run();
}