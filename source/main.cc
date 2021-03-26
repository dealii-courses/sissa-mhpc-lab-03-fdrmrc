#include "step-3.h"

int
main()
{
  deallog.depth_console(2);
  Step3 laplace_problem;
  laplace_problem.modify_bdary_cond = true; // add indicator 1
  laplace_problem.modify_bdary_data =
    true; // set value -1 to the bdary with indicator -1
  laplace_problem.source_term = 1.0;
  std::cout << laplace_problem.source_term << "\n";
  laplace_problem.run();

  return 0;
} 