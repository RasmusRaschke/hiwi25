#include "structures.hpp"
#include "systems.hpp"
#include "types.hpp"

CentralGravity::CentralGravity(double mu_) : mu(mu_) {}

vec3 CentralGravity::acceleration(const State& s, double) const {
    double r = s.r.norm();
    return -mu*s.r / (r*r*r);
}

double CentralGravity::potential(const State& s, double) const {
    return -mu / s.r.norm();
}