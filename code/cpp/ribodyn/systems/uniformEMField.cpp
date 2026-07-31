#include "systems.hpp"
#include "utility.hpp"

#include "systems.hpp"
#include "utility.hpp"

UniformEMField::UniformEMField(const vec3& electricField_, const vec3& magneticField_) : E0(electricField_), B0(magneticField_){}

double UniformEMField::scalarPotential(const State& state, double) const {
    return -E0.dot(state.r);
}

vec3 UniformEMField::vectorPotential(const State& state, double) const{
    return 0.5 * B0.cross(state.r);
}

vec3 UniformEMField::scalarPotentialGradient(const State&, double) const{
    return -E0;
}

mat3 UniformEMField::vectorPotentialJacobian(const State&, double) const{
    return 0.5 * util::hat(B0);
}

vec3 UniformEMField::vectorPotentialTimeDerivative(const State&, double) const{
    return vec3::Zero();
}

vec3 UniformEMField::electricField(const State&, double) const{
    return E0;
}

vec3 UniformEMField::magneticField(const State&, double) const{
    return B0;
}

mat3 UniformEMField::magneticFieldJacobian(const State&, double) const{
    return mat3::Zero();
}