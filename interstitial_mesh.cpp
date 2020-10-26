#include "interstitial_mesh.hpp"
#include <algorithm>
#include <cmath>

// remove Coordinate and change to Eigen::Vector3d everywhere
std::vector<Eigen::Vector3d> make_grid_points(int in_a, int in_b, int in_c, const casmutils::xtal::Lattice& lattice)
{
    // returns mxnxl grid of potential Li sites
    std::vector<Eigen::Vector3d> grid;
    Eigen::VectorXd Divisions_in_a = Eigen::VectorXd::LinSpaced(in_a, 0, 1);
    Eigen::VectorXd Divisions_in_b = Eigen::VectorXd::LinSpaced(in_b, 0, 1);
    Eigen::VectorXd Divisions_in_c = Eigen::VectorXd::LinSpaced(in_c, 0, 1);
    // triple for loop to add in all the coordinate divisions
    for (int i = 0; i < in_a; i++)
    {
        for (int j = 0; j < in_b; j++)
        {
            for (int k = 0; k < in_c; k++)
            {
                Eigen::Vector3d temp_coord;
                temp_coord << Divisions_in_a(i), Divisions_in_b(j), Divisions_in_c(k);
                Eigen::Vector3d temp_coord_cart = convert_to_cartesian(lattice, temp_coord);
                grid.emplace_back(temp_coord_cart);
            }
        }
    }
    // make sure this comes out as cartesian
    return grid;
}

// returns the coordinates of interstitial sites within the radius of a certain atom. The intent is to take in intersitial/atom coords in
// cartesian  coords, but a radius in cartesian coords, and then return the interstitial sites within the radius in cartesian coords. Does
std::vector<Eigen::Vector3d> find_interstitials_within_radius(std::vector<Eigen::Vector3d>& interstitial_coordinates,
                                                              const Eigen::Vector3d& sphere_origin,
                                                              double radius,
                                                              const casmutils::xtal::Lattice& lattice)
{

    // subtract sphere origin and then do distance algorithim (ie use peridic compare)
    std::vector<Eigen::Vector3d> interstitials_within_radius;
    for (const auto& interstitial_coord : interstitial_coordinates)
    {
        Eigen::Vector3d normalized_interstitial_coord = interstitial_coord - sphere_origin;
        VectorPeriodicCompare_f compare_interstitials_relative_to_single_site(Eigen::Vector3d(0, 0, 0), radius, lattice);
        if (compare_interstitials_relative_to_single_site(normalized_interstitial_coord) == true)
        {
            interstitials_within_radius.emplace_back(interstitial_coord);
        }
    }
    return interstitials_within_radius;
}

// Discard certain interstital grid points within certain radius
std::vector<Eigen::Vector3d> keep_reasonable_interstitial_gridpoints(const std::vector<Eigen::Vector3d>& total_interstitial_coordinates, const std::vector<Eigen::Vector3d>& interstitial_coordinates_to_discard, double tol, casmutils::xtal::Lattice& lattice)
{
	std::vector<Eigen::Vector3d> remaining_interstitials;
	for (const Eigen::Vector3d& interstitial_coord : total_interstitial_coordinates)
	{
		VectorPeriodicCompare_f compare_total_to_discarding_interstitials(interstitial_coord, tol, lattice);
		if (std::find_if(interstitial_coordinates_to_discard.begin(), interstitial_coordinates_to_discard.end(),
 compare_total_to_discarding_interstitials)==interstitial_coordinates_to_discard.end())
		{
			remaining_interstitials.emplace_back(interstitial_coord);
		}
	}
	return remaining_interstitials;
}


std::vector<Eigen::Vector3d> make_orbit(const Eigen::Vector3d& coordinate,
					const std::vector<casmutils::sym::CartOp>& symops,  
					const casmutils::xtal::Lattice& lattice)
{
	double tol=1e-5;
	std::vector<Eigen::Vector3d> orbit;
        for (const casmutils::sym::CartOp& symop: symops)
	{
            Eigen::Vector3d transformedcoord = symop * coordinate;
            VectorPeriodicCompare_f make_sure_no_repeated_coords(transformedcoord, tol, lattice);
	    if (std::find_if(orbit.begin(), orbit.end(), make_sure_no_repeated_coords)==orbit.end())
	//    //Need condition to check if the created coord is not exactly the same as another coord VectorPeriodicCompare_f maybe
			    {
			    orbit.emplace_back(transformedcoord);	
			    }    
        }
	return orbit;
}

// Indexes vectors of the different asymmetric orbits of the system. This cam further be broken into the asymmetric atoms in the system
// For reference these use the interstitial coordinates after the old ones are taken out
std::vector<int>
label_by_symmetrical_equivalence(const std::vector<Eigen::Vector3d>& coordinates,
                                  const std::vector<casmutils::sym::CartOp>& symops,
				  const casmutils::xtal::Lattice& lattice,
                                  double tol)
{
    int label=0;
    std::vector<int> coordinate_tags(coordinates.size(),-1);
    for(int ix = 0; ix<(coordinate_tags.size()); ++ix)
    {
        if(coordinate_tags[ix]!=-1)
        {
            continue;
        }

        const Eigen::Vector3d& coord=coordinates[ix];
        std::vector<Eigen::Vector3d> coord_orbit=make_orbit(coord, symops,lattice);

        for(int ixx=ix; ixx<coordinate_tags.size(); ++ixx)
        {
            if(coordinate_tags[ixx]!=-1)
            {
                continue;
            }

            VectorPeriodicCompare_f equals_ixx_coord(coordinates[ixx],tol,lattice);
	    if(std::find_if(coord_orbit.begin(),coord_orbit.end(),equals_ixx_coord)!=coord_orbit.end())
            {
                coordinate_tags[ixx]=label;	
            }
        }
	++label;
    }

    for(auto l : coordinate_tags)
    {
        assert(l!=-1);
    }
	
    return coordinate_tags;
}

std::vector<std::vector<Eigen::Vector3d>>
bin_by_symmetrical_equivalence(const std::vector<Eigen::Vector3d>& coordinates,
                                  const std::vector<casmutils::sym::CartOp>& symops,
		                  const casmutils::xtal::Lattice& lattice,
                                  double tol)
{
	std::vector<int> coordinate_tags= label_by_symmetrical_equivalence(coordinates, symops, lattice, tol);
	//goes through length of coordinate list and connects coordinate to orbit	
	int num_orbits=*std::max_element(coordinate_tags.begin(), coordinate_tags.end())+1;
	std::vector<std::vector<Eigen::Vector3d>> orbit_container;
	orbit_container.resize(num_orbits);
       	for (int i=0; i<coordinates.size(); i++)
	{
		Eigen::Vector3d temp_coord=coordinates[i];
		int label=coordinate_tags[i];
		orbit_container[label].emplace_back(temp_coord);
	}
	return orbit_container;	
}

//function for taking each of the points/"orbits" you have in the mesh, applying the factor group and getting the full orbit
std::vector<std::vector<Eigen::Vector3d>> apply_factor_group_on_each_bin(const std::vector<Eigen::Vector3d>& coordinates, const std::vector<casmutils::sym::CartOp>& symops, const casmutils::xtal::Lattice& lattice, double tol)
{
	std::vector<std::vector<Eigen::Vector3d>> partial_orbit_container=bin_by_symmetrical_equivalence(coordinates, symops, lattice, tol);
	std::vector<std::vector<Eigen::Vector3d>> full_orbit_container;
	for (const std::vector<Eigen::Vector3d>& partial_orbit :partial_orbit_container) 
	{
		std::vector<Eigen::Vector3d> orbits_to_append=make_orbit(partial_orbit[0], symops, lattice);
		full_orbit_container.push_back(orbits_to_append);
        }
	return full_orbit_container;

}

std::vector<Eigen::Vector3d> make_asymmetric_unit(const std::vector<Eigen::Vector3d>& complete_structure_basis,
                                                  const std::vector<casmutils::sym::CartOp>& symops,
						  const casmutils::xtal::Lattice& lattice,
                                                  double tol)
{
	std::vector<Eigen::Vector3d> asymmetric_unit_collated;
	std::vector<std::vector<Eigen::Vector3d>> total_orbits=bin_by_symmetrical_equivalence(complete_structure_basis, symops, lattice, tol);
	for (const auto& orbit: total_orbits)
	{
	if (orbit.size()>0)
	  {  
			asymmetric_unit_collated.push_back(orbit[0]);
	  }
	else
		std::cout<<"orbit of size zero here!"<<std::endl;	
	}
	return asymmetric_unit_collated;
}

