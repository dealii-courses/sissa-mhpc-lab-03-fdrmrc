#include "step-3.h"


void
convergence_mean_value();
int
main()
{
  deallog.depth_console(2);
  Step3 laplace_problem;
  laplace_problem.modify_bdary_cond = false; // add indicator 1
  laplace_problem.modify_bdary_data =
    false; // set value -1 to the bdary with indicator -1
  laplace_problem.source_term = 1.0;
  laplace_problem.l_shaped    = false;
  // std::cout << laplace_problem.source_term << "\n";
  laplace_problem.run();
  convergence_mean_value();
  return 0;
}



void
convergence_mean_value()
{
  std::cout << "\n"
            << "STARTING CONVERGENCE ANALYSIS FOR MEAN VALUE"
            << "\n";
  const unsigned int n_successive_refs = 7;

  for (unsigned int i = 1; i < n_successive_refs; ++i)
    {
      Step3 laplace_problem_conv;
      // square, dirichlet boundary conditions on eaach face
      laplace_problem_conv.modify_bdary_cond = false;
      laplace_problem_conv.modify_bdary_data = false;
      laplace_problem_conv.l_shaped          = false;
      laplace_problem_conv.n_global_refs     = i;

      laplace_problem_conv.run();
    }
}