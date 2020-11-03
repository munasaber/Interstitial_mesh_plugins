#include "../../CASMcode/include/casm/external/Eigen/Core"
#include "../../CASMcode/include/casm/external/Eigen/Dense"
#include <algorithm>
#include <casmutils/definitions.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/symmetry.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include "interstitial_mesh.hpp"
#include <filesystem>
#include <iostream>
#include <fstream>
#include <vector>
#include <CLI/CLI.hpp>
#include "CLI/App.hpp"
#include "CLI/Formatter.hpp"
#include "CLI/Config.hpp"


std::vector<Eigen::Vector3d> coordinates_to_remove(std::vector<double> radii, std::vector<Eigen::Vector3d> grid_dimensions, std::vector<std::vector<Eigen::Vector3d>> original_structure_eigen_coordinates, casmutils::xtal::Lattice lattice, double tol)
{
	std::vector<Eigen::Vector3d> coordinate_removal_list;
	for (int i=0; i<radii.size(); i++)
	{
		for (const Eigen::Vector3d& original_atom_site : original_structure_eigen_coordinates[i])
		{
			std::vector<Eigen::Vector3d> interstitials_to_remove_per_radius=find_interstitials_within_radius(grid_dimensions, original_atom_site, radii[i], lattice);	
			//do within radii for all atoms then compare within radii atoms to all atoms
			for (const Eigen::Vector3d& removal_atom: interstitials_to_remove_per_radius)
			{
				VectorPeriodicCompare_f test_coord(removal_atom, tol, lattice);
				if(std::find_if(grid_dimensions.begin(), grid_dimensions.end(), test_coord)!= grid_dimensions.end())
				{
					coordinate_removal_list.push_back(removal_atom);
			
				}
			}
		}
	}

	return coordinate_removal_list;
}


std::vector<std::vector<Eigen::Vector3d>> find_orbits_for_specific_mesh_dimensions(casmutils::fs::path structurepath, int int_a, int int_b, int int_c, std::vector<double> radii, std::vector<std::string> atomtypes,  double tol)
{
	casmutils::xtal::Structure original_structure=casmutils::xtal::Structure::from_poscar(structurepath);
	std::vector<casmutils::xtal::Site> structure_basis_sites=original_structure.basis_sites();
	casmutils::xtal::Lattice lattice=original_structure.lattice();
	std::vector<Eigen::Vector3d> grid_dimensions=make_grid_points(int_a, int_b, int_c, lattice);
	std::vector<std::vector<Eigen::Vector3d>> original_structure_eigen_coordinates;
	original_structure_eigen_coordinates.resize(atomtypes.size());
	for (int i=0; i<atomtypes.size(); i++)
	{
		for (const auto& site : structure_basis_sites)
		{
			if (site.label()==atomtypes[i])
			{
				original_structure_eigen_coordinates[i].push_back(site.cart());
			}
		}	
	}
        std::vector<Eigen::Vector3d> coordinate_removal_list=coordinates_to_remove(radii, grid_dimensions, original_structure_eigen_coordinates, lattice, tol);
	std::vector<Eigen::Vector3d> seived_grid_dimensions=keep_reasonable_interstitial_gridpoints(grid_dimensions, coordinate_removal_list, tol, lattice); 
	std::vector<casmutils::sym::CartOp> symops=make_factor_group(original_structure, tol);
	std::vector<std::vector<Eigen::Vector3d>> full_orbit_container=apply_factor_group_on_each_bin(seived_grid_dimensions, symops, lattice, tol);
        return full_orbit_container;
}	


//user inputted mesh function declaration

int main(int argc, char* argv[]) { 
	double tol=1E-4;
	CLI::App app{"This app takes in either a coordinate or the dimensions of a mesh and outputs equivalent orbits"};
	casmutils::fs::path structurepath;
	CLI::Option* structure_path=app.add_option("-s, --structure", structurepath, "Please input the base structure")-> required();
	
	std::vector<std::string> atomtypes;
	CLI::Option* atomtypes_path= app.add_option("-a, --atomtypes", atomtypes, "Please put the names of the different atoms in the structure") -> allow_extra_args(true);
	std::vector<double> distances;
	CLI::Option* distances_path= app.add_option("-d, --distances", distances, "Please input a list of the distances between the interstitial coordinates and each atom type in the base structure in Angstrom");
	std::vector<int> mesh;
	CLI::Option* mesh_path= app.add_option("-m, --mesh", mesh, "Please input the mesh dimensions in 'x,y,z' that you would like to examine in the base structure")->expected(3);	
	std::string outpath;
	CLI::Option* out_path= app.add_option("-o, --output", outpath, "Output path name");	
        structure_path->check(CLI::ExistingFile);
	CLI11_PARSE(app, argc, argv);
	std::cout<<"The chosen POSCAR is"<<structurepath<< std::endl;
	casmutils::xtal::Structure original_structure=casmutils::xtal::Structure::from_poscar(structurepath);

	std::cout<<"The chosen atoms are ";
	std::cout<<std::endl;
	for (const auto& my_atom: atomtypes)
	{
		std::cout<<my_atom<<std::endl;	
	}
	std::cout<<std::endl;
	std::cout<< "The chosen distances are ";
	std::cout<<std::endl;
        for (const auto& my_distance: distances)
	{
		std::cout<<my_distance<<std::endl;	
	}
	std::cout<<std::endl;
	std::cout<<"The chosen mesh is" << std::endl;	
	
	for (const auto& mesh_dimensions : mesh)
	{
		std::cout<<mesh_dimensions<<std::endl;	
	}
	std::cout << std::endl;
	
	std::vector<std::vector<Eigen::Vector3d>> all_orbits=find_orbits_for_specific_mesh_dimensions(structurepath, mesh[0], mesh[1], mesh[2], distances, atomtypes, tol); 
        int i=0;	
        std::ofstream my_outpath_file(outpath);
	for (const auto& orbit: all_orbits)
	{
		i++;
		my_outpath_file<<"These are the coordinates within orbit number "<<i<<std::endl;
			for (const auto& coordinate: orbit)
			{
				my_outpath_file<<coordinate.transpose()<<std::endl;
			}
			my_outpath_file<<std::endl;
	}
	my_outpath_file.close();
	return 0; 
}

