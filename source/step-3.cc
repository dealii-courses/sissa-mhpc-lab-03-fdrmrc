/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2020 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 *
 * Authors: Wolfgang Bangerth, 1999,
 *          Guido Kanschat, 2011
 *          Luca Heltai, 2021
 */

#include "step-3.h"

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

using namespace dealii;

Step3::Step3()
  : l_shaped{false}
  , modify_bdary_cond{false}
  , modify_bdary_data{false}
  , source_term{1.0} // default value for source term is 1
  , n_global_refs{5}
  , fe(1)
  , dof_handler(triangulation)
{}



void
Step3::make_grid()
{
  if (l_shaped)
    {
      GridGenerator::hyper_L(triangulation, -1.0, 1.0, true);
    }
  else
    {
      GridGenerator::hyper_cube(triangulation, -1.0, 1.0);
      if (modify_bdary_cond)
        {
          (triangulation.begin_active())
            ->face(0)
            ->set_boundary_id(1); // set a face with a bdary indicator = 1
        }
    }


  triangulation.refine_global(n_global_refs);
  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
}



void
Step3::setup_system()
{
  dof_handler.distribute_dofs(fe);
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit(sparsity_pattern);
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}



void
Step3::assemble_system()
{
  QGauss<2>          quadrature_formula(fe.degree + 1);
  FEValues<2>        fe_values(fe,
                        quadrature_formula,
                        update_values | update_gradients | update_JxW_values);
  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell_matrix = 0;
      cell_rhs    = 0;
      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          for (const unsigned int i : fe_values.dof_indices())
            for (const unsigned int j : fe_values.dof_indices())
              cell_matrix(i, j) +=
                (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                 fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                 fe_values.JxW(q_index));           // dx
          for (const unsigned int i : fe_values.dof_indices())
            cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                            source_term *                       // f(x_q)
                            fe_values.JxW(q_index));            // dx
        }
      cell->get_dof_indices(local_dof_indices);
      for (const unsigned int i : fe_values.dof_indices())
        for (const unsigned int j : fe_values.dof_indices())
          system_matrix.add(local_dof_indices[i],
                            local_dof_indices[j],
                            cell_matrix(i, j));
      for (const unsigned int i : fe_values.dof_indices())
        system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           Functions::ConstantFunction<2>(0.0),
                                           boundary_values);

  if (modify_bdary_data && modify_bdary_cond && l_shaped)
    {
      VectorTools::interpolate_boundary_values(
        dof_handler,
        1,
        Functions::ConstantFunction<2>(0.0),
        boundary_values); // used to test bdary indicators

      VectorTools::interpolate_boundary_values(
        dof_handler,
        2,
        Functions::ConstantFunction<2>(0.0),
        boundary_values); // used to test bdary indicators

      VectorTools::interpolate_boundary_values(
        dof_handler, 3, Functions::ConstantFunction<2>(0.0), boundary_values);

      VectorTools::interpolate_boundary_values(
        dof_handler, 4, Functions::ConstantFunction<2>(0.0), boundary_values);

      VectorTools::interpolate_boundary_values(
        dof_handler, 5, Functions::ConstantFunction<2>(0.0), boundary_values);
    }
  else if (modify_bdary_cond && modify_bdary_data)
    {
      // VectorTools::interpolate_boundary_values(
      // dof_handler, 0, Functions::ConstantFunction<2>(1.0), boundary_values);
    }
  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);
}



void
Step3::solve()
{
  SolverControl            solver_control(2000, 1e-12);
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
}



void
Step3::output_results() const
{
  DataOut<2> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();
  std::ofstream output("solution.vtk");
  data_out.write_vtk(output);

  std::cout << "Mean value: "
            << VectorTools::compute_mean_value(dof_handler,
                                               QGauss<2>(fe.degree + 1),
                                               solution,
                                               0)
            << std::endl;
}



void
Step3::run()
{
  make_grid();
  setup_system();
  assemble_system();
  solve();
  output_results();
}
