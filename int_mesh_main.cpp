#include "../../CASMcode/include/casm/external/Eigen/Core"
#include "../../CASMcode/include/casm/external/Eigen/Dense"
#include <algorithm>
#include <casmutils/definitions.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/symmetry.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include "interstitial_mesh.hpp"
#include <exception>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <vector>
#include <CLI/CLI.hpp>
#include "CLI/App.hpp"
#include "CLI/Formatter.hpp"
#include "CLI/Config.hpp"


std::vector<Eigen::Vector3d> coordinates_to_remove(const std::vector<double>& radii, std::vector<Eigen::Vector3d> grid_dimensions, const std::vector<std::vector<Eigen::Vector3d>>& original_structure_eigen_coordinates, const casmutils::xtal::Lattice& lattice, double tol)
{
	//remove sites based on being too close to an existing atom in the structure, based on each atoms cutoff radii
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

std::vector<Eigen::Vector3d> get_brought_within_coordinates(const std::vector<Eigen::Vector3d>& coordinate_vectors, const casmutils::xtal::Lattice& lattice)
{
	std::vector<Eigen::Vector3d> brought_within_coordinate_vectors;
	for (const auto& coordinate_vector: coordinate_vectors)
	{
		brought_within_coordinate_vectors.emplace_back(casmutils::xtal::cartesian_to_fractional(bring_within_lattice(coordinate_vector, lattice), lattice));
	}
	return brought_within_coordinate_vectors;
}	

std::vector<std::vector<Eigen::Vector3d>> find_orbits_for_specific_mesh_dimensions(casmutils::fs::path structurepath, int int_a, int int_b, int int_c, const std::vector<double>& radii, const std::vector<std::string>& atomtypes,  double tol)
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
	//std::vector<std::vector<Eigen::Vector3d>> full_orbit_container=apply_factor_group_on_each_bin(get_brought_within_seived_grid_dimensions(seived_grid_dimensions, lattice), symops, lattice, tol);
        return full_orbit_container;
}	


std::vector<Eigen::Vector3d> UniqueVectors(const std::vector<Eigen::Vector3d>& grid_dimensions, const casmutils::xtal::Lattice& lattice, double tol)
{
	std::vector<Eigen::Vector3d> seived_dimensions;
	for (const auto& point: grid_dimensions)
	{
		VectorPeriodicCompare_f test_coord(point, tol, lattice);
		if (std::find_if(seived_dimensions.begin(), seived_dimensions.end(), test_coord)==seived_dimensions.end())
		{
			seived_dimensions.emplace_back(point);
		}	
	}
	return seived_dimensions;
}


//Function for outputting a poscar that includes all sites in the poscar
casmutils::xtal::Structure final_mesh_structure(casmutils::fs::path& structurepath, std::string interstitialtype, int int_a, int int_b, int int_c, const std::vector<double>& radii, const std::vector<std::string>& atomtypes, double tol)
{
	//get structure and lattice
	casmutils::xtal::Structure original_structure= casmutils::xtal::Structure::from_poscar(structurepath);
	casmutils::xtal::Lattice lattice= original_structure.lattice();	
	
	//organize atomtypes in structure to match radii cutoffs to atomtypes
	std::vector<casmutils::xtal::Site> structure_basis_sites=original_structure.basis_sites();
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

	//construct the list of coordinates to be removed for being within the cutoff radius for each specietype
	std::vector<Eigen::Vector3d> coordinate_removal_list=coordinates_to_remove(radii, grid_dimensions, original_structure_eigen_coordinates, lattice, tol);
	std::vector<Eigen::Vector3d> seived_grid_points=keep_reasonable_interstitial_gridpoints(grid_dimensions, coordinate_removal_list, tol, lattice);
        std::vector<Eigen::Vector3d> unique_grid_points=UniqueVectors(seived_grid_points, lattice, tol);
       	for (const auto& coordinate : unique_grid_points)
	{
		structure_basis_sites.emplace_back(bring_within_lattice(coordinate, lattice), interstitialtype);	
	}
	casmutils::xtal::Structure output_structure(lattice, structure_basis_sites);
	return output_structure;
}


void print_dilute_poscar(casmutils::fs::path& structurepath, std::string interstitialtype, int int_a, int int_b, int int_c, const std::vector<double>& radii, const std::vector<std::string>& atomtypes, double tol, casmutils::fs::path outpath_file)
{
	//get structure and lattice
	casmutils::xtal::Structure original_structure= casmutils::xtal::Structure::from_poscar(structurepath);
	casmutils::xtal::Lattice lattice= original_structure.lattice();	
	
	//organize atomtypes in structure to match radii cutoffs to atomtypes
	std::vector<casmutils::xtal::Site> original_basis_sites=original_structure.basis_sites();
	std::vector<Eigen::Vector3d> grid_dimensions=make_grid_points(int_a, int_b, int_c, lattice);
	std::vector<std::vector<Eigen::Vector3d>> original_structure_eigen_coordinates;
	original_structure_eigen_coordinates.resize(atomtypes.size());
	for (int i=0; i<atomtypes.size(); i++)
	{
		for (const auto& site : original_basis_sites)
		{
			if (site.label()==atomtypes[i])
			{
				original_structure_eigen_coordinates[i].push_back(site.cart());
			}
		}	
	}

	//construct the list of coordinates to be removed for being within the cutoff radius for each specietype
	std::vector<Eigen::Vector3d> coordinate_removal_list=coordinates_to_remove(radii, grid_dimensions, original_structure_eigen_coordinates, lattice, tol);
	std::vector<Eigen::Vector3d> seived_grid_points=keep_reasonable_interstitial_gridpoints(grid_dimensions, coordinate_removal_list, tol, lattice);
        std::vector<Eigen::Vector3d> unique_grid_points=UniqueVectors(seived_grid_points, lattice, tol);
	int i=0;
	for (const auto& point: unique_grid_points)
	{
		std::vector<casmutils::xtal::Site> structure_basis_sites=original_structure.basis_sites();
		structure_basis_sites.emplace_back(bring_within_lattice(point, lattice), interstitialtype);
		casmutils::xtal::Structure output_structure(lattice, structure_basis_sites);
		std::string name=(outpath_file.string()+"_"+std::to_string(i));
		casmutils::fs::path enumerated_output_file=name;
		casmutils::xtal::write_poscar(output_structure, enumerated_output_file);
		i++;
	}

}


//user inputted mesh function declaration
int main(int argc, char* argv[]) { 
	double tol=1E-4;
	CLI::App app{"This app takes in either a coordinate or the dimensions of a mesh and outputs equivalent orbits"};
	casmutils::fs::path structurepath;
	CLI::Option* structure_path=app.add_option("-s, --structure", structurepath, "Please input the base structure")-> required();
	
	std::vector<std::string> atomtypes;
	CLI::Option* atomtypes_path= app.add_option("-a, --atomtypes", atomtypes, "Please put the names of the different atoms in the structure") -> allow_extra_args(true) -> required();
	std::vector<double> distances;
	CLI::Option* distances_path= app.add_option("-d, --distances", distances, "Please input a list of the distances between the interstitial coordinates and each atom type in the base structure in Angstrom") -> required();
	std::vector<int> mesh;
	CLI::Option* mesh_path= app.add_option("-m, --mesh", mesh, "Please input the mesh dimensions in 'x,y,z' that you would like to examine in the base structure")->expected(3) -> required();	
	std::string outpath;
	CLI::Option* out_path= app.add_option("-o, --output", outpath, "Output path name for orbits in the structure");
	std::string interstitialtype;
	CLI::Option* interstitial_type = app.add_option("-i, --interstitial", interstitialtype, "Name of interstitial atom type");
	std::string outputposcar;
	CLI::Option* output_poscar = app.add_option("-p, --output_poscar", outputposcar, "Output POSCAR with all interstitial atoms"); 	
	std::string diluteposcar;
	CLI::Option* dilute_poscar=app.add_option("--dilute_poscar", diluteposcar, "Output each poscar corresponding to a possible grid point");
	
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
	if (* out_path)
	{
		std::vector<std::vector<Eigen::Vector3d>> all_orbits=find_orbits_for_specific_mesh_dimensions(structurepath, mesh[0], mesh[1], mesh[2], distances, atomtypes, tol); 
        	int i=0;	
        	std::ofstream my_outpath_file(outpath);
		for (const auto& orbit: all_orbits)
		{
			i++;
			std::vector<Eigen::Vector3d> brought_within_orbit=get_brought_within_coordinates(orbit, original_structure.lattice());
			my_outpath_file<<"These are the coordinates within orbit number "<<i<<std::endl;
			for (const auto& coordinate: brought_within_orbit)
			{
				my_outpath_file<<coordinate.transpose()<<std::endl;
			}
			my_outpath_file<<std::endl;
		}
		my_outpath_file.close();
	}
	std::cout<<"After the outpath file loop"<<std::endl;
	if (* output_poscar)
	{
		std::cout<<"In the output poscar loop";
		//prints a poscar that includes all coordinates from the mesh, set them to the inierstitial site type and ouput the combination as a poscar
	        casmutils::xtal::write_poscar(final_mesh_structure(structurepath, interstitialtype, mesh[0], mesh[1], mesh[2], distances, atomtypes, tol), outputposcar);
	}
	if (* dilute_poscar)
	{
		std::cout<<"Printing a poscar for each valid grid point";
		print_dilute_poscar(structurepath, interstitialtype, mesh[0], mesh[1], mesh[2], distances, atomtypes, tol, diluteposcar);
	}
	return 0; 
}

