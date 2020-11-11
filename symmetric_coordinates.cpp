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
void print_out_poscar_from_user_specified_coordinated(Eigen::Vector3d cartesian_coordinate, casmutils::fs::path structurepath, std::string interstitialtype, double tol, casmutils::fs::path outpath)
{
	auto original_structure=casmutils::xtal::Structure::from_poscar(structurepath);
	std::vector<casmutils::xtal::Site> structure_basis_sites=original_structure.basis_sites();
	casmutils::xtal::Lattice lattice=original_structure.lattice();
	std::vector<casmutils::sym::CartOp> factor_group=make_factor_group(original_structure, tol);
	std::vector<Eigen::Vector3d> orbit_from_user_input_coordinate=make_orbit(cartesian_coordinate, factor_group, lattice); 
        for (const auto& coordinate_vector: orbit_from_user_input_coordinate)
	{
	 	const casmutils::xtal::Coordinate coordinate=casmutils::xtal::Coordinate(coordinate_vector);
		structure_basis_sites.emplace_back((coordinate).bring_within(lattice), interstitialtype);	
	}
	casmutils::xtal::Structure output_structure(lattice, structure_basis_sites);
	casmutils::xtal::write_poscar(output_structure, outpath);	
}



//user inputted mesh function declaration

int main(int argc, char* argv[]) { 
	double tol=1E-4;
	CLI::App app{"This app takes in either a coordinate or the dimensions of a mesh and outputs equivalent orbits"};
	casmutils::fs::path structurepath;
	CLI::Option* structure_path=app.add_option("-s, --structure", structurepath, "Please input the base structure")-> required();
	std::vector<double> Cartesian_coordinate;
	CLI::Option* Cartesian_coordinate_path= app.add_option("-c, --Cartesian_coordinate", Cartesian_coordinate, "Please input a sinlge point that you would like to find the orbit of within the structure in Cartesian coordinates");
	std::vector<double> Fractional_coordinate;
	CLI::Option* Fractional_coordinate_path= app.add_option("-f, --Fractional_coordinate", Fractional_coordinate, "Please input a sinlge point that you would like to find the orbit of within the structure in Fractional coordinates");
	std::string interstitialtype;
	CLI::Option* interstitialtype_path= app.add_option("-i, --interstitialtype", interstitialtype, "Please put the names of the different atoms in the structure")->required();
	
	casmutils::fs::path outpath;
	CLI::Option* out_path= app.add_option("-o, --output", outpath, "Output path name");	
        structure_path->check(CLI::ExistingFile);
	out_path->check(CLI::ExistingFile);
	CLI11_PARSE(app, argc, argv);
	std::cout<<"The chosen POSCAR is"<<structurepath<< std::endl;
	casmutils::xtal::Structure original_structure=casmutils::xtal::Structure::from_poscar(structurepath);
        Eigen::Vector3d vector_coordinate;
	if (Cartesian_coordinate.size()==3)
	{
		vector_coordinate<<Cartesian_coordinate[0], Cartesian_coordinate[1], Cartesian_coordinate[2];	
	        print_out_poscar_from_user_specified_coordinated(vector_coordinate, structurepath, interstitialtype, tol, outpath);
		
	}
	else if (Fractional_coordinate.size()==3)
	{
		casmutils::xtal::Lattice lattice=original_structure.lattice();
	        vector_coordinate<<Fractional_coordinate[0], Fractional_coordinate[1], Fractional_coordinate[2];
		print_out_poscar_from_user_specified_coordinated(convert_to_cartesian(lattice, vector_coordinate), structurepath, interstitialtype, tol, outpath);
	}
	else 
	{
		std::cout<<"The input of the coordinate is not in a standard form (Vector of size 3)"<<std::endl;
	}
	return 0; 
}

