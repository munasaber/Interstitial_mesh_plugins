#ifndef COORDINATE_H
#define COORDINATE_H
#include "../../CASMcode-dev/include/casm/external/Eigen/Core"
#include "../../CASMcode-dev/include/casm/external/Eigen/Dense"
#include <casmutils/xtal/lattice.hpp>
#include <string>
#include <vector>
#include <casmutils/xtal/symmetry.hpp>


//Compares two vectors without taking into account the periodicity of the lattice
//Vectors are meant to indicate the positions of an atomic site
struct VectorCompare_f
{
	VectorCompare_f(const Eigen::Vector3d& vector, double prec);
	bool operator()(const Eigen::Vector3d& other) const;
	private:
	Eigen::Vector3d m_vector;
	double m_precision;
};

//added from lattice.hpp from original creation of the vectorfunctions file. 
//Consider moving it to the lattice function
Eigen::Vector3d convert_to_cartesian(const casmutils::xtal::Lattice& lattice, const Eigen::Vector3d& frac_coord); 

//added from lattice.hpp from original creation of the vectorfunctions file. 
//Consider moving it to the lattice function
Eigen::Vector3d convert_to_fractional(const casmutils::xtal::Lattice& lattice, const Eigen::Vector3d& cart_coord); 

//Compares two vectors, taking into account the periodicity of the lattice
//Vectors are meant to indicate the positions of an atomic site
struct VectorPeriodicCompare_f
{
    VectorPeriodicCompare_f(const Eigen::Vector3d& vector, double prec, const casmutils::xtal::Lattice& unit_cell);
    bool operator()(const Eigen::Vector3d& other) const;
private:
    Eigen::Vector3d m_vector;
    double m_precision;
    casmutils::xtal::Lattice m_lattice;
};


//Operator for the mutliplication of a symmetry operation and a vector 
Eigen::Vector3d operator*(const casmutils::sym::CartOp& transformation, const Eigen::Vector3d& vector);

#endif
