#ifndef INTERSTITIAL_MESH_H
#define INTERSTITIAL_MESH_H

#include "casmutils/xtal/lattice.hpp"
#include "casmutils/xtal/site.hpp"
#include "casmutils/xtal/symmetry.hpp"
#include "casmutils/xtal/structure.hpp"
#include <vector>
#include "vectorfunctions.hpp"
//Next tools we need:
//Geometric center of mass for cluster
//Find asymetric unit of a basis*
//Find all sites within a radius*


//makes grid points of an interstitial atoms of the specified dimensions using a certain lattice
std::vector<Eigen::Vector3d> make_grid_points(int in_a, int in_b, int in_c, const casmutils::xtal::Lattice& lattice);

//This function finds which interstitial atoms are within a certain radius of a particular atom. 
std::vector<Eigen::Vector3d> find_interstitials_within_radius(std::vector<Eigen::Vector3d>& interstitial_coordinates, const Eigen::Vector3d& sphere_origin, double radius, const casmutils::xtal::Lattice& lattice);

//This is used in conjunction with the make_gridpoints function. 
//It prevents overlaps of the interstitial atoms with the atoms in the base structure
std::vector<Eigen::Vector3d> keep_reasonable_interstitial_gridpoints(const std::vector<Eigen::Vector3d>& total_interstitial_coordinates, const std::vector<Eigen::Vector3d>& interstitial_coordinates_to_discard, double tol, casmutils::xtal::Lattice& lattice);

//This makes an orbit (or a set of symmetrically equivalent atoms)
std::vector<Eigen::Vector3d> make_orbit(const Eigen::Vector3d& coordinates, const std::vector<casmutils::sym::CartOp>& symops, const casmutils::xtal::Lattice& lattice);


//This function labels atoms by their respective orbits
std::vector<int>
label_by_symmetrical_equivalence(const std::vector<Eigen::Vector3d>& coordinates,
                                  const std::vector<casmutils::sym::CartOp>& symops,
                                  const casmutils::xtal::Lattice& lattice,
                                  double tol);

//This function is used to bin sets of atoms into different orbits
std::vector<std::vector<Eigen::Vector3d>>
bin_by_symmetrical_equivalence(const std::vector<Eigen::Vector3d>& coordinates,
                                  const std::vector<casmutils::sym::CartOp>& symops,
                                  const casmutils::xtal::Lattice& lattice,
                                  double tol);

//Each bin of orbis has it's factor group applied 
std::vector<std::vector<Eigen::Vector3d>> apply_factor_group_on_each_bin(const std::vector<Eigen::Vector3d>& coordinates, const std::vector<casmutils::sym::CartOp>& symops, const casmutils::xtal::Lattice& lattice, double tol);

//finds the asymmetric unit of a ser of atoms through application of the factor group
std::vector<Eigen::Vector3d> make_asymmetric_unit(const std::vector<Eigen::Vector3d>& complete_structure_basis, 
                                  const std::vector<casmutils::sym::CartOp>& symops,
				  const casmutils::xtal::Lattice& lattice, 
			          double tol);
#endif 
