#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>



#define EARTH_GRAV 9.80665;

using namespace std;
using vec3 = Eigen::Vector3d;
using quat = Eigen::Quaterniond;
using mat3 = Eigen::Matrix3d;

struct Parameters {
    double M;       //Mass 
    double R;       //Radius
    double g;       //Gravitational Acceleration
    double phi;     //Incline Angle, tan phi = alpha
    double I;       //inertia scalar 
    double nu;      //friction coefficient
    double eta;     //air resistance coefficient
    vec3 n;         //normal vector field
    vec3 mu_body;   //magnetic moment in body frame
    vec3 B_inert;    //magnetic field in inertial frame
};

struct State {
    vec3 r;
    vec3 Omega;
    quat q;
};

template<typename>
struct always_false : std::false_type {};

std::string remove_comment(std::string line) {
    size_t pos = line.find("#");
    if (pos != std::string::npos) {
        line = line.substr(0, pos);
    }
    line.erase(0, line.find_first_not_of(" \t"));
    line.erase(line.find_last_not_of(" \t") + 1);
    return line;
}

template<typename T>
T get_input_variable(const std::string& line) {
    std::string cleaned_line = remove_comment(line);
    try {
        if constexpr (std::is_same_v<T, double>) {
            return std::stod(cleaned_line);
        } else if constexpr (std::is_same_v<T, int>) {
            return std::stoi(cleaned_line);
        } else {
            static_assert(always_false<T>::value, "Unsupported type");
        }
    } catch (const std::exception& e) {
        std::cerr << "Error parsing input variable: " << e.what() << "\n";
        exit(1);
    }
}

template<typename T>
T get_next_input_variable(){
    std::string line;
    if (!std::getline(std::cin, line)) {
        std::cerr << "Error reading input variable: End of input reached\n";
        exit(1);
    }
    return get_input_variable<T>(line);
}

static quat exponential(const vec3 &rotvec){
    double theta = rotvec.norm();
    if (theta < 1e-15) return quat::Identity();
    vec3 axis = rotvec / theta;
    double a = cos(0.5 * theta);
    double b = sin(0.5 * theta);
    return quat(a, b * axis.x(), b * axis.y(), b * axis.z());
}

static vec3 rotate(const quat &q, const vec3 &v){
    return q * v;
}

static vec3 torque(const quat &q, const vec3 &mu, const vec3 &B){
    vec3 mu_init = rotate(q, mu);
    return mu_init.cross(B);
}

static vec3 getOmega(const vec3 &Omega, const quat &q, const Parameters &p){
    const vec3 e3(0.,0.,1.);
    double N = p.M * p.g * p.n.dot(e3);
    vec3 tang = -(e3 - (e3.dot(p.n)) * p.n);
    double tnorm = tang.norm();
    if (tnorm > 1e-15) tang /= tnorm;
    else tang.setZero();
    double vpar = (p.R * Omega.cross(p.n)).dot(tang);
    double sgn = 0.0;
    if (vpar > 1e-15) sgn = 1.0;
    else if (vpar < -1e-15) sgn = -1.0;
    vec3 Fr = p.nu * N * sgn * tang;
    vec3 tau = torque(q, p.mu_body, p.B_inert) - p.R * (p.n.cross(Fr));
    vec3 dotOmega_t = - ((5 * p.g) / (7 * p.R)) * p.n.cross(e3) + (5 / (7 * p.R * p.R * p.M)) * (tau - (p.n.dot(tau)) * p.n);
    double dotOmega_n = (p.n.dot(tau)) / p.I;
    return (dotOmega_t + dotOmega_n * p.n)- p.eta * Omega;
}

static vec3 getr(const vec3 &Omega, const Parameters &p){
    return p.R * Omega.cross(p.n);
}


static State integrator(const State &s, double dt, const Parameters &p){
    vec3 Omega1 = s.Omega;
    //RKMK Steps
    quat q1 = s.q;
    vec3 k1_O = getOmega(Omega1, q1, p);
    vec3 k1_r = getr(Omega1, p);
    vec3 Omega2 = s.Omega + 0.5 * dt * k1_O;
    quat q2 = exponential(0.5 * dt * Omega1) * s.q;
    vec3 k2_O = getOmega(Omega2, q2, p);
    vec3 k2_r = getr(Omega2, p);
    vec3 Omega3 = s.Omega + 0.5 * dt * k2_O;
    quat q3 = exponential(0.5 * dt * Omega2) * s.q;
    vec3 k3_O = getOmega(Omega3, q3, p);
    vec3 k3_r = getr(Omega3, p);
    vec3 Omega4 = s.Omega + dt * k3_O;
    quat q4 = exponential(dt * Omega3) * s.q;
    vec3 k4_O = getOmega(Omega4, q4, p);
    vec3 k4_r = getr(Omega4, p);  
    //RKMK Evaluation
    vec3 Omega_update = s.Omega + (dt / 6.0) * (k1_O + 2.0 * k2_O + 2.0 * k3_O + k4_O);
    vec3 r_update = s.r + (dt / 6.0) * (k1_r + 2.0 * k2_r + 2.0 * k3_r + k4_r);
    //Quaternion Update
    quat dq = exponential(dt * 0.5 * (s.Omega + Omega_update));
    quat q_update = dq * s.q;
    return State{r_update, Omega_update, q_update};
}

int main(){
    double t = 0.0;
    double dt = get_next_input_variable<double>();
    double t_end = get_next_input_variable<double>();

    Parameters p;
    p.M = get_next_input_variable<double>();
    p.R = get_next_input_variable<double>();
    p.g = get_next_input_variable<double>();
    p.phi = get_next_input_variable<double>() * M_PI / 180.0;
    p.n = vec3(0.0, -sin(p.phi), cos(p.phi));
    double mu_x = get_next_input_variable<double>();
    double mu_y = get_next_input_variable<double>();
    double mu_z = get_next_input_variable<double>();
    p.mu_body = vec3(mu_x, mu_y, mu_z);
    double B_x = get_next_input_variable<double>();
    double B_y = get_next_input_variable<double>();
    double B_z = get_next_input_variable<double>();
    p.B_inert = vec3(B_x, B_y, B_z);
    double start_spin = get_next_input_variable<double>() * M_PI / 180.0; // convert to radians per second
    //p.mu_body = vec3(0.0, 0.0, 0.0); 
    //p.B_inert = vec3(0.0, 0.0, 0.0);
    p.nu = get_next_input_variable<double>();
    p.eta = get_next_input_variable<double>();
    p.I = (2.0 / 5.0) * p.M * p.R * p.R;
    std::cout << p.mu_body << "\n";
    std::cout << p.B_inert << "\n";

    State s;
    s.r = vec3(0.0, 0.0, p.R*(1.0 - p.n.dot(vec3(0,0,1))));
    s.Omega = p.n * start_spin;
    std::cout << p.n << "\n";
    s.q = quat::Identity();

    int steps = int(t_end / dt);
    
    ofstream out("data.csv");
    out << scientific << setprecision(9);
    out << "t,x,y,z,Ox,Oy,Oz,q0,q1,q2,q3,";
    out << "vx, vy, vz, T_trans, T_rot, U_gr, U_em, E, q_norm\n";
    for (int i=0; i <= steps; ++i){
        if ((i % 10) == 0){
            vec3 v = getr(s.Omega, p);
            double T_trans = 0.5 * p.M * v.squaredNorm();
            double T_rot = 0.5 * p.I * s.Omega.squaredNorm();
            double U_gr = p.M * p.g * s.r.z();
            double U_em = - (rotate(s.q, p.mu_body)).dot(p.B_inert);
            double E = T_trans + T_rot + U_gr + U_em;
            double q_norm = s.q.norm();
            out << t << "," << s.r.x() << "," << s.r.y() << "," << s.r.z() << ","
            << s.Omega.x() << "," << s.Omega.y() << "," << s.Omega.z() << ","
            << s.q.w() << "," << s.q.x() << "," << s.q.y() << "," << s.q.z() << ","
            << v.x() << "," << v.y() << "," << v.z() << "," << T_trans << "," << T_rot << "," << U_gr << "," << U_em << "," << E << "," << q_norm << "\n";
        }
        s = integrator(s, dt, p);
        t += dt;
    }
    out.close();
    cout << "Done";
    return 0;
}
