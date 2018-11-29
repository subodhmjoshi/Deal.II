/*

	This program gives matrices (Mass and Stiffness) for tensor-product cells. 
	These values can be used for validating the code that you have written for computing these matrices.
	To run,
	1. Go to the main() function (at the end of the file), and enter values for order etc.
	2. Save and exit
	3. compile and run like the usual deal.ii file
	
 */


#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <fstream>
#include <iostream>
#include <vector>

using namespace dealii;



template<int dim>
class DGMat
{
public:
  DGMat (int order, int intFlag, int refinementLevels);

  void run ();

private:
  int order;
  int K; //order+1 for inexact integration, order+2 for exact integration
  int RL; // refinementLevels
  void make_grid ();
  void setup_system ();
  void assemble_system ();

  Triangulation<dim>     triangulation;
  FE_Q<dim>              fe;
  DoFHandler<dim>        dof_handler;

};


template<int dim>
DGMat<dim>::DGMat (int order, int intFlag, int refinementLevels)
  :
  fe (order), 
  dof_handler (triangulation)
{
	assert(intFlag == 0 or intFlag == 1);
	//0: inexact
	//1: exact
	this->order = order;
	this->K = order + 1 + intFlag;
	this->RL = refinementLevels;
}


template<int dim>
void DGMat<dim>::make_grid ()
{
  std::cout << "Domain boundary [-1,1] in each direction. Make suitable changes in DGMat<dim>::make_grid()\n";
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (RL);

  std::cout << "Order of the polynomial: "
            << this->order
            << std::endl;
}

template<int dim>
void DGMat<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);
  std::cout << "\nNumber of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;
}


template<int dim>
void DGMat<dim>::assemble_system (){
	// LGL points and roots
	QGaussLobatto<dim>  quadrature_formula(this->K);

	//LG points and roots
	//QGauss<dim>  quadrature_formula(this->K);

	FEValues<dim> fe_values (fe, quadrature_formula,
			       update_values | update_gradients | update_JxW_values | update_quadrature_points);

	const unsigned int   dofs_per_cell = fe.dofs_per_cell;
	const unsigned int   n_q_points    = quadrature_formula.size();


	std::cout << "\nIntegration Type: (0: inexact, 1:exact):  " << this->K - (this->order+1) << std::endl;
	std::cout << "\nNumber of quadrature points: " << n_q_points << std::endl;


	std::vector<std::vector < double > > M (dofs_per_cell, std::vector<double> (dofs_per_cell)); //Mass matrix
	std::vector<std::vector < double > > D (dofs_per_cell, std::vector<double> (dofs_per_cell)); //Differentiation matrix (weak form)


	typename std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
	typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
	typename DoFHandler<dim>::active_cell_iterator endc = dof_handler.end();

	int cellCounter = 0;
	for (; cell!=endc; ++cell){
		std::cout << "\n\n\nCell number: " << cellCounter << " (note: counting started from 0)\n";
		cellCounter ++;
		
		fe_values.reinit (cell);

	      

        	// Computing Mass matrix
	      for (unsigned int i=0; i<dofs_per_cell; ++i){
		for (unsigned int j=0; j<dofs_per_cell; ++j){
		  for (unsigned int q_index=0; q_index<n_q_points; ++q_index){
			  M[i][j] += (fe_values.shape_value (i, q_index) *
					   fe_values.shape_value (j, q_index) *
					   fe_values.JxW (q_index));
		  };
		};
	      };

	      
	      //printing mass matrix
	      std::cout << "Mass matrix: \n";
	      for (unsigned int i=0; i<dofs_per_cell; i++){
		      for (unsigned int j=0; j<dofs_per_cell ; j++){
			      std::cout << M[i][j] << " ";
			      //resetting the Mass matrix for the next cell computation:
			      M[i][j] = 0.0;
		      };
		      std::cout << std::endl;
	      };

	      // Mass matrix has to be rearranged for typical applications.
	      // The reason: In deal.ii, the polynomials are arranged as first, last, second, third ... 
	      // i.e. the DOF locations are such that the endpoints are arranged first (for Gauss lobatto points) then the rest
	      // example: For 1D element, the 2nd order interpolation will have DOF locations as: (-1,1,0)
	      // Thus, the second row and the second column have to be shifted to the last row and the last column respectively.

	      std::cout << "\nIMPORTANT: Kindly take note that all the cell matrices have the rows and columns mixed up. This is because the DOF locations are [-1, 1, rest] instead of [-1,rest,1] in 1D. \n" << std :: endl;


	      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
        	// Computing Diff in x direction matrix
	      for (unsigned int i=0; i<dofs_per_cell; ++i){
		for (unsigned int j=0; j<dofs_per_cell; ++j){
		  for (unsigned int q_index=0; q_index<n_q_points; ++q_index){
			  D[i][j] += (fe_values.shape_grad (i, q_index)[0] *
					   fe_values.shape_value (j, q_index) *
					   fe_values.JxW (q_index));
		  };
		};
	      };

	      
	      std::cout << "Dx matrix: \n";
	      for (unsigned int i=0; i<dofs_per_cell; i++){
		      for (unsigned int j=0; j<dofs_per_cell ; j++){
			      std::cout << D[i][j] << " ";
			      //Resetting value of D for next computations
			      D[i][j] = 0.0;
		      };
		      std::cout << std::endl;
	      };
	      std::cout << "\nIMPORTANT: Kindly take note that all the cell matrices have the rows and columns mixed up. This is because the DOF locations are [-1, 1, rest] instead of [-1,rest,1] in 1D. \n" << std :: endl;


	     if (dim > 1){ 
		      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			// Computing Diff in y direction matrix
		      for (unsigned int i=0; i<dofs_per_cell; ++i){
			for (unsigned int j=0; j<dofs_per_cell; ++j){
			  for (unsigned int q_index=0; q_index<n_q_points; ++q_index){
				  D[i][j] += (fe_values.shape_grad (i, q_index)[1] *
						   fe_values.shape_value (j, q_index) *
						   fe_values.JxW (q_index));
			  };
			};
		      };

		      std::cout << "Dy matrix: \n";
		      for (unsigned int i=0; i<dofs_per_cell; i++){
			      for (unsigned int j=0; j<dofs_per_cell ; j++){
				      std::cout << D[i][j] << " ";
				      //Resetting value of D for next computations
				      D[i][j] = 0.0;
			      };
			      std::cout << std::endl;
		      };
		      std::cout << "\nIMPORTANT: Kindly take note that all the cell matrices have the rows and columns mixed up. This is because the DOF locations are [-1, 1, rest] instead of [-1,rest,1] in 1D. \n" << std :: endl;



		      if (dim >2){
			      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
				// Computing Diff in z direction matrix
			      for (unsigned int i=0; i<dofs_per_cell; ++i){
				for (unsigned int j=0; j<dofs_per_cell; ++j){
				  for (unsigned int q_index=0; q_index<n_q_points; ++q_index){
					  D[i][j] += (fe_values.shape_grad (i, q_index)[2] *
							   fe_values.shape_value (j, q_index) *
							   fe_values.JxW (q_index));
				  };
				};
			      };

			      
			      std::cout << "Dz matrix: \n";
			      for (unsigned int i=0; i<dofs_per_cell; i++){
				      for (unsigned int j=0; j<dofs_per_cell ; j++){
					      std::cout << D[i][j] << " ";
					      //Resetting value of D for next computations
					      D[i][j] = 0.0;
				      };
				      std::cout << std::endl;
			      };
			      std::cout << "\nIMPORTANT: Kindly take note that all the cell matrices have the rows and columns mixed up. This is because the DOF locations are [-1, 1, rest] instead of [-1,rest,1] in 1D. \n" << std :: endl;
		      };
	     };
	};
}





template<int dim>
void DGMat<dim>::run ()
{
  make_grid ();
  setup_system ();
  assemble_system ();
}


int main ()
{
  deallog.depth_console (2);

  int intFlag = 0; //0: inexact;  1: exact
  int order = 1; // polynomial order
  int refinementLevels = 0; //number of levels of refinement of the hypbercube
 
  
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
  std::cout << "1D Matrices:\n";
  const int dim1 = 1; // number of spatial dimensions

  DGMat <dim1> Matrices1(order, intFlag, refinementLevels);
  Matrices1.run ();

  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
  std::cout << "\n\n\n2D Matrices:\n";
  const int dim2 = 2; // number of spatial dimensions

  DGMat <dim2> Matrices2(order, intFlag, refinementLevels);
  Matrices2.run ();

  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
  std::cout << "\n\n\n3D Matrices:\n";
  const int dim3 = 3; // number of spatial dimensions

  DGMat <dim3> Matrices3(order, intFlag, refinementLevels);
  Matrices3.run ();

  return 0;
}
