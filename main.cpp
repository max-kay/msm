#include "builder.cpp"
#include <iostream>
#include <vector>
#define EPSILON_0 (double)8.8541878188e-12
#define MY_0 (double)1.25663706127e-6
using namespace std;

// corrFactorVelocity
class Sim {
  public:
    SimulationParameters params;
    vector<v3> particle_e_field;
    vector<v3> particle_e_dipol;
    vector<v3> particle_h_field;
    vector<v3> particle_pos;
    vector<v3> particle_velocity;
    vector<v3> particle_direction;
    vector<v3> particle_direction_velocity;

    Sim(SimulationParameters params) : params(params) {
        // todo
    }

    double find_delta_t() {
        double max_value = 0;
        for (auto elem : particle_velocity) {
            double norm_velocity = elem.get_length();

            if (norm_velocity > max_value) {
                max_value = norm_velocity;
            }
        }
        double delta_t =
            params.corrFactorVelocity * params.dragRadius / max_value;

        return delta_t;
    };
    void update_position(double delta_t) {
        for (int i = 0; i < params.numberOfParticles; i++) {
            particle_pos[i] =
                (particle_pos[i] + particle_velocity[i] * delta_t) %
                params.lengthSimulationCube;
        }
    }

    void update_direction(double delta_t) {
        for (int i = 0; i < params.numberOfParticles; i++) {
            particle_direction[i] =
                particle_direction[i] +
                delta_t *
                    particle_direction_velocity[i]; // periodic BC are not
                                                    // required here right?
        }
    }

    void update_pos_velocity() {
        double R_radius =
            cbrt(pow(params.longSemiaxesAB, 2) * params.chiEffShortAxisC);
        double h_prefactor =
            3 * MY_0 *
            pow(params.magMomentDensityParticle * params.volumeParticle, 2) /
            (4 + M_PI);
        double e_prefactor =
            3 / (EPSILON_0 * params.relPermittivityMatrix * 2 * M_PI);
        double r_prefactor =
            3 * MY_0 *
            pow(params.magMomentDensityParticle * params.volumeParticle, 2) /
            (2 * M_PI * pow(2 * R_radius, 4));
        for (int i = 0; i < params.numberOfParticles; i++) {
            v3 this_h_force = v3(0, 0, 0);
            v3 this_e_force = v3(0, 0, 0);
            v3 this_r_force = v3(0, 0, 0);
            for (int j = 0; i < params.numberOfParticles; i++) {
                if (j == i) {
                    continue;
                }
                v3 r_ji = (particle_pos[i] - particle_pos[j]) %
                          params.lengthSimulationCube; // distance with PBC
                v3 r_ji_hat = r_ji.get_direction();
                this_h_force = // chekc the beackets
                    this_h_force +
                    h_prefactor * (1 / pow(r_ji.get_length(), 4)) *
                        ((((particle_direction[j].dot(particle_direction[i])) -
                           5 * (r_ji_hat.dot(particle_direction[j])) *
                               (r_ji_hat.dot(particle_direction[i]))) *
                          r_ji_hat) +
                         ((r_ji_hat.dot(particle_direction[i])) *
                          particle_direction[j]) +
                         (r_ji_hat.dot(particle_direction[j]) *
                          particle_direction[i]));
                this_e_force =
                    this_e_force +
                    e_prefactor / pow(r_ji.get_length(), 4) *
                        ((((particle_e_dipol[j].dot(particle_e_dipol[i])) -
                           5 * (r_ji_hat.dot(particle_e_dipol[j])) *
                               (r_ji_hat.dot(particle_e_dipol[i]))) *
                          r_ji_hat) +
                         (((r_ji_hat.dot(particle_e_dipol[i])) *
                           particle_e_dipol[j]) +
                          (r_ji_hat.dot(particle_e_dipol[j])) *
                              particle_e_dipol[i]));
                this_r_force =
                    this_r_force +
                    r_prefactor *
                        exp(-params.corrFactorRepulsiveForce *
                            ((r_ji.get_length() / 2 * R_radius) - 1)) *
                        r_ji_hat;
            } // second for
            particle_velocity[i] =
                1 / params.dragCoeffTransl *
                (this_e_force + this_h_force +
                 this_r_force); // not adding to it is correcft right?
        } // first for
    } // entire func

    void update_direc_velocity() {
        // todo
    }

    void update_e_field() {
        double prefactor = 1 / (4 * M_PI * params.relPermittivityMatrix *
                                EPSILON_0); // NO DIELECTRIC CONST. epslion_0

        for (int i = 0; i < params.numberOfParticles; i++) {
            // these are here to calculate the distance from our i-th
            // particle to the j-th particle
            v3 this_e_field = v3(0, 0, 0);
            for (int j = 0; j < params.numberOfParticles; j++) {
                if (j == i) {
                    continue;
                } // dont devide by 0
                v3 r_ji = (particle_pos[i] - particle_pos[j]) %
                          params.lengthSimulationCube; // distance with PBC
                v3 r_ji_hat = r_ji.get_direction();
                this_e_field = // calculate the electric field of the particle
                    this_e_field +
                    (prefactor / (pow(r_ji.get_length(), 3))) *
                        (3.0 * (r_ji_hat *
                                (particle_direction[j].dot(r_ji_hat))) -
                         particle_direction[j]);
            }
            particle_e_field[i] = this_e_field; // update the list
        }
    }
    void update_h_field() { // check math pls otherwise analogous to
                            // update_e_field function
        for (int i = 0; i < params.numberOfParticles; i++) {
            v3 this_h_field = v3(0, 0, 0);
            double prefactor = params.magMomentDensityParticle *
                               params.volumeParticle / (4 * M_PI);
            for (int j = 0; j < params.numberOfParticles; i++) {
                if (j == i) {
                    continue;
                }
                v3 r_ji_h = (particle_pos[i] - particle_pos[j]) %
                            params.lengthSimulationCube;
                v3 r_ji_h_hat = r_ji_h.get_direction();
                this_h_field =
                    this_h_field +
                    prefactor * (1 / pow(r_ji_h.get_length(), 3)) *
                        (3 * (particle_direction[j].dot(r_ji_h_hat)) *
                             r_ji_h_hat -
                         particle_direction[j]);
            }
            particle_h_field[i] = this_h_field;
        }
    }
    void update_el_dipole() { // pls check
        double prefactor = params.volumeParticle * EPSILON_0;
        double chi_diff = params.chiEffShortAxisC - params.chiEffLongAxesAB;
        for (int i; i < params.numberOfParticles; i++) {
            v3 left_term = params.chiEffLongAxesAB * particle_e_field[i];
            v3 right_term = chi_diff *
                            (particle_direction[i].dot(particle_e_field[i])) *
                            particle_direction[i];
            particle_e_dipol[i] = prefactor * (left_term + right_term);
        }
    }
};

int main() {
    std::cout << "Hello World" << std::endl;
    return 0;
}
