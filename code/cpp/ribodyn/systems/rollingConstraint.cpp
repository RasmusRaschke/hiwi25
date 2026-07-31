#include "structures.hpp"
#include "systems.hpp"
#include "utility.hpp"

ConstraintData RollingConstraint::evaluate(const State& s, double) const {
    ConstraintData c;
    c.A.resize(3,6);
    const mat3 R = s.q.toRotationMatrix();
    const vec3 nBody = R.transpose() * normal;
    c.A.leftCols<3>() = mat3::Identity();
    c.A.rightCols<3>() = radius * R * util::hat(nBody);
    c.b = Eigen::VectorXd::Zero(3);
    c.gamma = Eigen::VectorXd::Zero(3);
    return c;
}