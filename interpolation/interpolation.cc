#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>

using namespace dealii;

enum Type {GL, GLob, Gch};

class parameter{
	public:
		parameter ();
		double xl,xu;
		int n, order;
		Type type;
		void print();
};

parameter::parameter(){}

void parameter::print(){
	std::cout << "\n______________________________________________________________" << std::endl;
	switch (type) {
		case 0:
			std::cout << "Gauss-Legendre Quadrature Formula"<< std::endl;
			break;
		case 1:
			std::cout << "Gauss-Lobatto Quadrature Formula"<< std::endl;
			break;
		case 2:
			std::cout << "Gauss-Chebichev Quadrature Formula"<< std::endl;
			break;
	};
	std::cout << "Xl: " << xl << "  Xu: " << xu << std::endl;
	std::cout << "n: " << n << "   order: " << order << std::endl;
	std::cout << "______________________________________________________________" << std::endl;
}




class Interpolation{
	public:
		Interpolation(parameter param);
		void interpolate();
		void print();
	private:
		void build_grid();
		void initialize();
		void compute_norm();
		void plot_error();

		Triangulation<1>  triangulation;
		FE_Q<1>           FE;
		DoFHandler<1>     DOF;

		Vector <double>   Q_exact;
		Vector <double>   Q_num;
		Vector <double>   Q_sample;
		Vector <double>   X;


		Type type;
		double xl,xu;
		unsigned int n, order;



};

Interpolation::Interpolation(parameter param):
	FE(param.order), //Here FE knows the degree of the polynomial. FE contains information about finite element space.
	DOF(triangulation), //Here DOF object is associated with the triangulation object.
	xl(param.xl),xu(param.xu),order(param.order),type(param.type),n(param.n)
{}

void Interpolation::build_grid(){
	GridGenerator::hyper_cube (triangulation,xl,xu); //triangulation gets the geometrical and physical structure
	std::cout << "\nNumber of active cells: "
		  << triangulation.n_active_cells()
		  << std::endl;
	DOF.distribute_dofs(FE); //Depending on the FE object (remember that the FE already knows the DOF of the problem) DOF distributed;
	//note: FE links triangulation and DOF objects together.
	std::cout << "\nNumber of Sampling points: "
		  << DOF.n_dofs()
		  << std::endl;

	std::cout << "\nNumber of Plotting points: "
		  << n
		  << std::endl;

	Q_exact.reinit(n);        //solution plotted at n points
	Q_num.reinit(n);          //solution plotted at n points
	X.reinit(n);              //x coordinates at n points
	Q_sample.reinit(order+1); //Sampling points in the Lagrange formula
}


void Interpolation::initialize(){

	//here exact solution (and later the numerical solution) has been found at Quadrature points
	QGauss<1>    Q_exact_points(n);
	//(order+1) polynomials will be computed for n quadrature points
	// however, at this points, we are not interested in the created polynomials
	// we just want to get the quadrature points' geometrical location to set the
	// initial conditions. 
	// So only thing that will be used is the quadrature_points() from this formula.
	FEValues<1>  fe_values_exact(FE,Q_exact_points,update_values|update_quadrature_points);

	const unsigned int ndof = FE.dofs_per_cell;
	const unsigned int qsize = Q_exact_points.size();
	DoFHandler<1>::active_cell_iterator 
		cell = DOF.begin_active(), 
		endc = DOF.end();


	// setting initial solution on all plotting points
	for (;cell<endc;++cell){
		fe_values_exact.reinit(cell); //map the canonical element to the physical cell
		for (int i=0; i<qsize;++i){
			double x = fe_values_exact.quadrature_point(i)[0]; //x value of the quad. point
			Q_exact(i) = std::sin(M_PI * x);
			Q_num(i) = 0;
			X(i) = x;
		};
	};

	// setting value of Q on all sampling points (these are y_i corresponding to x_i in lagrange formula)
	// Note: L_2 projection method can be used for this purpose, however, here `direct' approach is used.
	// First step is to get the location of the sampling points.
	// This should depend on the DOF distribution over the cell. Following functions extract that 
	// information
	MappingQ1<1> mapping; //define mapping
	std::vector<Point<1> > support_points(DOF.n_dofs());  //create array for storing points
	DoFTools::map_dofs_to_support_points(mapping,DOF,support_points); //map points 

	// now directly assign the values based on the x location 
	for (int i=0;i<ndof;i++){
		double x =  support_points[i][0] ; //x coordinate only
		Q_sample[i] = std::sin(M_PI * x);
	};
}
	

void Interpolation::interpolate(){
	build_grid();
	initialize();

	// now we need to create the interpolation polynomials 
	// we also create the identical quad. point. as the exact solution locations.
	// thus we compute the lagrange polynomials at the quad. points. based on the sampling points
	QGauss<1>    Interpolation_points(n);
	FEValues<1>  fe_values(FE,Interpolation_points,update_values|update_quadrature_points|update_JxW_values);

	const unsigned int ndof = FE.dofs_per_cell;
	const unsigned int qsize = Interpolation_points.size();
	DoFHandler<1>::active_cell_iterator 
		cell = DOF.begin_active(), 
		endc = DOF.end();

	// standard deal (deal.ii !!)
	for (;cell<endc;++cell){
		fe_values.reinit(cell);
		for (int q=0; q<qsize; ++q){
			// now we are on a quadrature point
			for (int i=0; i<ndof; i++){
				// ith polynomial is considered
				Q_num(q) += Q_sample(i)*fe_values.shape_value(i,q);  //shape value is the ith Lagrange polynomial computed at the qth quadrature point.
			};
		};
	};
}

	

void Interpolation::print(){
	std::cout << "\n\n" << std::endl;
	std::cout << "x \t\t Exact \t\t Numerical" << std::endl;
	std::cout << "------------------------------------------------------" << std::endl;
	for (int i=0; i<n; i++){
		std::cout << X(i) << " \t" << Q_exact(i) << " \t" << Q_num(i) << std::endl;
	}
}


int main(){
	parameter p1;
	p1.xl = 0;
	p1.xu = 1;
	p1.order = 2; //note: polynomial order = number of datapoints - 1
	p1.n = 25;
	p1.type = GL;

	Interpolation I(p1);

	I.interpolate();

	I.print();

	return 0;
}




// Todo:
// Error computation based on the exact and numerical solutions
// Implementation of the Chebichev and Lobatto points
// Automated Error plotting routine
