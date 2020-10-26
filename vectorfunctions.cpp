#include "vectorfunctions.hpp"

// Vector Compare functor
//Compares two vectors without taking into account the periodicity of the lattice
//Vectors are meant to indicate the positions of an atomic site
VectorCompare_f::VectorCompare_f(const Eigen::Vector3d& vector, double prec) : m_vector(vector), m_precision(prec) {}
bool VectorCompare_f::operator()(const Eigen::Vector3d& other) const
{
    if (m_vector.isApprox(other, m_precision))
    {
            return true;
    }
    return false;
}


//Converts vector signifying site coordinates to fractional coordinates 
Eigen::Vector3d convert_to_fractional(const casmutils::xtal::Lattice& lattice, const Eigen::Vector3d& cart_coord) 
{
   Eigen::Vector3d frac_coord= lattice.column_vector_matrix().inverse()*cart_coord;
   return frac_coord;
}



//Converts vector signifying site coordinates to cartesian coordinates 
Eigen::Vector3d convert_to_cartesian(const casmutils::xtal::Lattice& lattice, const Eigen::Vector3d& frac_coord) 
{
   Eigen::Vector3d cart_coord= lattice.column_vector_matrix()*frac_coord;
   return cart_coord;
}


//Compares two vectors, taking into account the periodicity of the lattice
//Vectors are meant to indicate the positions of an atomic site
VectorPeriodicCompare_f::VectorPeriodicCompare_f(const Eigen::Vector3d& vector, double prec, const casmutils::xtal::Lattice& unit_cell) : m_vector(vector), m_precision(prec), m_lattice(unit_cell){
} 
bool VectorPeriodicCompare_f::operator()(const Eigen::Vector3d& other) const 
{	
    Eigen::Vector3d vector1 = convert_to_fractional(m_lattice, m_vector);
    Eigen::Vector3d vector2 = convert_to_fractional(m_lattice, other);
    Eigen::Vector3d distance_vector=vector1-vector2;
    
    for (int i=0; i<distance_vector.size();i++){
        distance_vector(i)=distance_vector(i)-std::round(distance_vector(i));
    }
    Eigen::Vector3d cartesian_distance_vector = convert_to_cartesian(m_lattice, distance_vector);
    return std::abs(cartesian_distance_vector.norm())<m_precision;
}


//Operator for the mutliplication of a symmetry operation and a vector 
Eigen::Vector3d operator*(const casmutils::sym::CartOp& transformation, const Eigen::Vector3d& vector)
{
    Eigen::Matrix3d transformed_matrix= (transformation.matrix);
    Eigen::Vector3d transformed = (transformed_matrix)*vector+transformation.translation;
    return transformed;
}

