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
#include <vector>
#include <CLI/CLI.hpp>
#include "CLI/App.hpp"
#include "CLI/Formatter.hpp"
#include "CLI/Config.hpp"

//user inputted coordinate function declaration
//prints out a poscar with the user specified coordinate and symmetrically equivalent atoms
//void print_out_poscar_from_user_specified_coordinated(Eigen::Vector3d coordinate, casmutils::fs::path structurepath, std::string interstitialtype, double tol, casmutils::fs::path outpath)
//{
//	auto original_structure=casmutils::xtal::Structure::from_poscar(structurepath);
//	std::vector<casmutils::xtal::Site> structure_basis_sites=original_structure.basis_sites();
//	casmutils::xtal::Lattice lattice=original_structure.lattice();
//	std::vector<casmutils::sym::CartOp> factor_group=make_factor_group(original_structure, tol);
//	std::vector<Eigen::Vector3d> orbit_from_user_input_coordinate=make_orbit(coordinate, factor_group, lattice); 
//        for (const auto& coordinate_vector: orbit_from_user_input_coordinate)
//	{
//	 	structure_basis_sites.emplace_back(casmutils::xtal::Coordinate(coordinate_vector), interstitialtype);	
//	}
//	casmutils::xtal::Structure output_structure(lattice, structure_basis_sites);
//	casmutils::xtal::write_poscar(output_structure, outpath);	
//}

void print_orbits_for_specific_mesh_dimensions(casmutils::fs::path structurepath, int int_a, int int_b, int int_c, std::vector<double> radius, std::vector<std::string> atomtype, std::string interstitialtype, double tol, casmutils::fs::path outpath)
{
	casmutils::xtal::Structure original_structure=casmutils::xtal::Structure::from_poscar(structurepath);
	std::vector<casmutils::xtal::Site> structure_basis_sites=original_structure.basis_sites();
	casmutils::xtal::Lattice lattice=original_structure.lattice();
	std::vector<Eigen::Vector3d> grid_dimensions=make_grid_points(int_a, int_b, int_c, lattice);
	std::vector<Eigen::Vector3d> coordinate_removal_list;
	std::vector<std::vector<Eigen::Vector3d>> original_structure_eigen_coordinates;
	for (int i=0; i<atomtype.size(); i++)
	{
		for (const auto& site : structure_basis_sites)
		{
			if (site.label()==atomtype[i])
				original_structure_eigen_coordinates[i].emplace_back(site.cart());
		}	
	}
	for (int i=0; i<radius.size(); i++)
	{
		for (const Eigen::Vector3d& original_atom_site : original_structure_eigen_coordinates[i])
		{
			std::vector<Eigen::Vector3d> interstitials_to_remove_per_radius=find_interstitials_within_radius(grid_dimensions, original_atom_site, radius[i], lattice);	
			//do within radius for all atoms then compare within radius atoms to all atoms
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
	std::vector<Eigen::Vector3d> seived_grid_dimensions=keep_reasonable_interstitial_gridpoints(grid_dimensions, coordinate_removal_list, tol, lattice); 
	std::vector<casmutils::sym::CartOp> symops=make_factor_group(original_structure, tol);
	std::vector<std::vector<Eigen::Vector3d>> full_orbit_container=apply_factor_group_on_each_bin(seived_grid_dimensions, symops, lattice, tol);

}	


//user inputted mesh function declaration

int main(int argc, char* argv[]) { 
	double tol=1E-4;
	CLI::App app{"This app takes in either a coordinate or the dimensions of a mesh and outputs equivalent orbits"};
	casmutils::fs::path structurepath;
	CLI::Option* structure_path=app.add_option("-s, --structure", structurepath, "Please input the base structure")-> required();
	
	//std::vector<std::string> atomtype;
	//CLI::Option* atomtype_path= app.add_option("-a, --atomtype", atomtype, "Please put the names of the different atoms in the structure") -> allow_extra_args(true);
	//std::vector<double> distances;
	//CLI::Option* distances_path= app.add_option("-d, --distances", distances, "Please input a list of the distances between the interstitial coordinates and each atom type in the base structure in Angstrom");
	//std::vector<double> mesh(3);
	//CLI::Option* mesh_path= app.add_option("-m, --mesh", mesh, "Please input the mesh dimensions in 'x,y,z' that you would like to examine in the base structure");	
	casmutils::fs::path outpath;
	CLI::Option* out_path= app.add_option("-o, --output", outpath, "Output path name");	
        structure_path->check(CLI::ExistingFile);
	CLI11_PARSE(app, argc, argv);
	std::cout<<"The chosen POSCAR is"<<structurepath<< std::endl;
	casmutils::xtal::Structure original_structure=casmutils::xtal::Structure::from_poscar(structurepath);

	//std::cout<<"The chosen atoms are ";
	//std::cout<<std::endl;
	//for (const auto& my_atom: atomtype)
	//{
	//	std::cout<<my_atom<<std::endl;	
	//}
	//std::cout<<std::endl;
	//std::cout<< "The chosen distances are ";
	//std::cout<<std::endl;
        //for (const auto& my_distance: distances)
	//{
	//	std::cout<<my_distance<<std::endl;	
	//}
	//std::cout<<std::endl;
	//std::cout<<"The chosen mesh is" << std::endl;	
	//
	//for (const auto& mesh_dimensions : mesh)
	//{
	//	std::cout<<mesh_dimensions<<std::endl;	
	//}
	return 0; 
}

