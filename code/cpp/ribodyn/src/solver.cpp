#include "solver.hpp"
#include "structures.hpp"
#include <stdexcept>

Solver::Solver(const MechanicalSystem& system, SolverMode mode) : system(system), mode(mode){};

mat6 Solver::massMatrix() const {
    mat6 M = mat6::Zero();
    M.topLeftCorner<3,3>() = system.body.mass * mat3::Identity();
    M.bottomRightCorner<3,3>() = system.body.inertia;
    return M;
};

ConstraintData Solver::assembleConstraints(const State& s,double t) const{
    std::vector<ConstraintData> cache;
    cache.reserve(system.constraints.size());
    int totalRows = 0;
    for (const auto& constraint : system.constraints){
        ConstraintData data = constraint->evaluate(s, t);
        const int rows = static_cast<int>(data.A.rows());
        if(data.A.cols() != 6){
            throw std::runtime_error(
                "Constraint matrix must have six columns."
            );
        }
        if(data.b.size() != rows){
            throw std::runtime_error(
                "Constraint b vector has incorrect size."
            );
        }
        if(data.gamma.size() != rows){
            throw std::runtime_error(
                "Constraint gamma vector has incorrect size."
            );
        }
        totalRows += rows;
        cache.push_back(std::move(data));
    }
    ConstraintData result;
    result.A.resize(totalRows, 6);
    result.b.resize(totalRows);
    result.gamma.resize(totalRows);
    int currentRow = 0;
    for(const ConstraintData& data : cache){
        const int rows = static_cast<int>(data.A.rows());
        result.A.block(currentRow, 0, rows, 6) = data.A;
        result.b.segment(currentRow, rows) = data.b;
        result.gamma.segment(currentRow, rows) = data.gamma;
        currentRow += rows;
    }
    return result;
}

Wrench Solver::assembleGeneralizedForce(const State& s, double t) const {
    Wrench total;
    for(const auto& force : system.forces)
        total += force->evaluate(system.body, s, t);
    for(const auto& field : system.gravityFields){
        total.force += system.body.mass * field -> acceleration(s,t);
    }
    for(const auto& field : system.emFields){
        const vec3 E = field->electricField(s,t);
        const vec3 B = field->magneticField(s,t);
        total.force += system.body.charge * (E + s.v.cross(B));
        const mat3 R = s.q.toRotationMatrix();
        const vec3 muInertial = R * system.body.magneticMoment;
        const mat3 JB = field->magneticFieldJacobian(s, t);
        total.force += JB.transpose() * muInertial;
        const vec3 BBody = R.transpose() * B;
        total.torque += system.body.magneticMoment.cross(BBody);
    }
    return total;
};

Wrench Solver::assembleLagrangianWrench(const State& s, double t) const {
    Wrench total;
    for(const auto& p : system.potentials){
        const Wrench grad = p->gradient(system.body, s, t);
        total.force -= grad.force;
        total.torque -= grad.torque;
    }
    for(const auto& field : system.gravityFields){
        total.force += system.body.mass * field->acceleration(s, t);
    }
    for(const auto& field : system.emFields){
        const vec3 gradPhi = field->scalarPotentialGradient(s, t);
        const mat3 JA = field->vectorPotentialJacobian(s, t);
        const vec3 dAdt = field->vectorPotentialTimeDerivative(s, t);
        total.force += system.body.charge * ((JA.transpose() - JA) * s.v - dAdt - gradPhi);
        const vec3 B = field->magneticField(s, t);
        const mat3 JB = field->magneticFieldJacobian(s, t);
        const mat3 R = s.q.toRotationMatrix();
        const vec3 muInertial = R * system.body.magneticMoment;
        total.force += JB.transpose() * muInertial;
        const vec3 BBody = R.transpose() * B;
        total.torque += system.body.magneticMoment.cross(BBody);
    }
    return total;
};

StateDerivative Solver::rhs(const State& s, double t){
    const Wrench wrench = mode == SolverMode::dAlembert ? assembleGeneralizedForce(s, t) : assembleLagrangianWrench(s, t);
    const ConstraintData constraints = assembleConstraints(s, t);
    const mat6 M = massMatrix();
    vec6 Q;
    Q.head<3>() = wrench.force;
    Q.tail<3>() = wrench.torque - s.Omega.cross(system.body.inertia * s.Omega);
    const int m = static_cast<int>(constraints.A.rows());
    vec6 acceleration;
    if (m == 0){
        lastLambda.resize(0);
        Eigen::LDLT<mat6> decomposition(M);
        if(decomposition.info() != Eigen::Success){
            throw std::runtime_error(
                "Mass matrix factorization failed."
            );
        }
        acceleration = decomposition.solve(Q);
    }
    else{
        Eigen::MatrixXd K = Eigen::MatrixXd::Zero(6 + m, 6 + m);
        K.topLeftCorner(6, 6) = M;
        K.topRightCorner(6, m) = -constraints.A.transpose();
        K.bottomLeftCorner(m, 6) = constraints.A;
        Eigen::VectorXd kktRhs(6 + m);
        kktRhs.head<6>() = Q;
        kktRhs.tail(m) = constraints.gamma;
        Eigen::FullPivLU<Eigen::MatrixXd>
            decomposition(K);
        if(decomposition.rank() < 6 + m){
            throw std::runtime_error(
                "Singular KKT system: constraints may be dependent."
            );
        }
        const Eigen::VectorXd solution =
            decomposition.solve(kktRhs);
        acceleration = solution.head<6>();
        lastLambda = solution.tail(m);
    }
    StateDerivative derivative;
    derivative.rDot = s.v;
    derivative.vDot = acceleration.head<3>();
    derivative.Omega = s.Omega;
    derivative.OmegaDot = acceleration.tail<3>();
    return derivative;
}

void Solver::projectVelocity(State& s, double t) const{
    const ConstraintData constraints = assembleConstraints(s, t);
    const int m = static_cast<int>(constraints.A.rows());
    if (m == 0)
        return;
    vec6 velocity;
    velocity.head<3>() = s.v;
    velocity.tail<3>() = s.Omega;
    const Eigen::VectorXd residual = constraints.b - constraints.A * velocity;
    if (residual.norm() < 1e-13)
        return;
    const mat6 M = massMatrix();
    Eigen::LDLT<mat6> massDecomposition(M);
    if(massDecomposition.info() != Eigen::Success){
        throw std::runtime_error(
            "Mass matrix factorization failed during projection."
        );
    }
    const Eigen::MatrixXd inverseMassAT = massDecomposition.solve(constraints.A.transpose());
    const Eigen::MatrixXd projectionMatrix = constraints.A * inverseMassAT;
    Eigen::FullPivLU<Eigen::MatrixXd>
        projectionDecomposition(projectionMatrix);
    if(projectionDecomposition.rank() < m){
        throw std::runtime_error(
            "Velocity projection failed: dependent constraints."
        );
    }
    velocity += inverseMassAT * projectionDecomposition.solve(residual);
    s.v = velocity.head<3>();
    s.Omega = velocity.tail<3>();
}

double Solver::constraintResidual(const State& s, double t) const{
    const ConstraintData constraints = assembleConstraints(s, t);
    if (constraints.A.rows() == 0)
        return 0.0;
    vec6 velocity;
    velocity.head<3>() = s.v;
    velocity.tail<3>() = s.Omega;
    return (constraints.A * velocity - constraints.b).norm();
}

const Eigen::VectorXd&
Solver::getLambda() const {
    return lastLambda;
}

Diagnostics Solver::diagnostics(const State& s, double t) const{
    Diagnostics d;
    const mat3 R = s.q.toRotationMatrix();
    d.T_trans = 0.5 * system.body.mass * s.v.squaredNorm();
    d.T_rot = 0.5 * s.Omega.dot(system.body.inertia * s.Omega);
    for(const auto& potential : system.potentials){
        d.U_generic += potential->value(system.body, s, t);
    }
    for(const auto& field : system.gravityFields){
        d.U_gr += system.body.mass * field->potential(s, t);
    }
    d.mu_body = system.body.magneticMoment;
    d.mu_world = R * d.mu_body;
    for(const auto& field : system.emFields){
        d.E_world += field->electricField(s, t);
        d.B_world += field->magneticField(s, t);
        d.U_em += system.body.charge * field->scalarPotential(s, t);
    }
    d.U_em -= d.mu_world.dot(d.B_world);
    d.E_body = R.transpose() * d.E_world;
    d.B_body = R.transpose() * d.B_world;
    d.E_total = d.T_trans + d.T_rot + d.U_generic + d.U_gr + d.U_em;
    d.constraintResidual = constraintResidual(s, t);
    d.quaternionNorm = s.q.norm();
    return d;
}