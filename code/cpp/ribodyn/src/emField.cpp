#include "structures.hpp"
#include "types.hpp"
#include "utility.hpp"
#include <cmath>

ElectromagneticField::ElectromagneticField(double spatialStep_, double timeStep_) : spatialStep(spatialStep_), timeStep(timeStep_) {}

vec3 ElectromagneticField::scalarPotentialGradient(const State& state, double t) const {
    vec3 gradient;
    for(int j=0; j<3; ++j){
        const double h = spatialStep * (1.0 + std::abs(state.r(j)));
        State plus = state;
        State minus = state;
        plus.r(j) += h;
        minus.r(j) -= h;
        gradient(j) = (scalarPotential(plus, t) - scalarPotential(minus, t)) / (2.0 * h);
    }
    return gradient;
}

mat3 ElectromagneticField::vectorPotentialJacobian(const State& state, double t) const {
    mat3 jacobian;
    for (int j=0; j<3; ++j){
        const double h = spatialStep * (1.0 + std::abs(state.r(j)));
        State plus = state;
        State minus = state;
        plus.r(j) += h;
        minus.r(j) -= h;
        jacobian.col(j) = (vectorPotential(plus, t) - vectorPotential(minus, t)) / (2.0 * h);
    }
    return jacobian;
}

vec3 ElectromagneticField::vectorPotentialTimeDerivative(const State& state, double t) const {
    const double h = timeStep * (1.0 + std::abs(t));
    return (vectorPotential(state, t+h) - vectorPotential(state, t-h)) / (2.0 * h);
}

vec3 ElectromagneticField::electricField(const State& state, double t) const {
    return -scalarPotentialGradient(state, t) - vectorPotentialTimeDerivative(state, t);
}

vec3 ElectromagneticField::magneticField(const State& state, double t) const {
    const mat3 jacobian = vectorPotentialJacobian(state, t);
    return util::vee(jacobian - jacobian.transpose());
}

mat3 ElectromagneticField::magneticFieldJacobian(const State& state, double t) const {
    mat3 jacobian;
    for(int j=0; j<3; ++j) {
        const double h = spatialStep * (1.0 + std::abs(state.r(j)));
        State plus = state;
        State minus = state;
        plus.r(j) += h;
        minus.r(j) -= h;
        jacobian.col(j) = (magneticField(plus, t) - magneticField(minus, t)) / (2.0 * h);
    }
    return jacobian;
}