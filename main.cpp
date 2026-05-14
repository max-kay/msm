#include "builder.cpp"
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>
#include <vector>

#define EL_CONVERGE_ITERATIONS 4
#define EPSILON_0 (double)8.8541878188e-12
#define MU_0 (double)1.25663706127e-6

using namespace std;

class Sim {
  public:
    mt19937 rng;
    SimulationParameters params;
    vector<v3> particle_e_field;
    vector<v3> particle_e_dipol;
    vector<v3> particle_h_field;
    vector<v3> particle_pos;
    vector<v3> particle_velocity;
    vector<v3> particle_direction;
    vector<v3> particle_direction_velocity;

    Sim(SimulationParameters params) : params(params) { //check order pls and if everything is here
        generate_positions();
        generate_direction();
        generate_e_field();
        generate_h_field();
        generate_direc_vel();
        generate_vel();
        first_update_e_dipole();

        // TODO:
        // Generate random position in the simulation cube, make sure they dont
        // overlap using eqRadius. -- DONE
        // Generate random unit vectors for the direction. -- DONE
        // Set the fields to their default values. -- H is done E needs updates
        // from aurel
        // update the electric dipole once fill both velocities with
        // 0 vectors -- DONE
    }
    void generate_positions() {
        srand(time(nullptr)); // seed the generator
        for (int i = 0; i < params.numberOfParticles; i++) {
            double rng_1 = (double)rand() / RAND_MAX - 0.5;
            double rng_2 = (double)rand() / RAND_MAX - 0.5;
            double rng_3 = (double)rand() / RAND_MAX - 0.5;

            particle_pos.push_back(v3(rng_1 * params.lengthSimulationCube,
                                      rng_2 * params.lengthSimulationCube,
                                      rng_3 * params.lengthSimulationCube));
            for (int j = 0; j < i; j++) {
                if (j == i) {
                    continue;
                }
                double r_ji = (particle_pos[j] - particle_pos[i]).get_length();
                if (r_ji < params.eqRadius) {
                    double rng_1 = (double)rand() / RAND_MAX - 0.5;
                    double rng_2 = (double)rand() / RAND_MAX - 0.5;
                    double rng_3 = (double)rand() / RAND_MAX - 0.5;
                    particle_pos[i] = v3(rng_1 * params.lengthSimulationCube,
                                         rng_2 * params.lengthSimulationCube,
                                         rng_3 * params.lengthSimulationCube);
                    j = -1; // reset the for loop bc we need to check all the
                            // values again if we change them
                }
            }
        }
    }

    void generate_direction() { // pls check if you find same math just to be safe
        for (int i = 0; i < params.numberOfParticles; i++) {
            double z =
                2.0 * ((double)rand() / RAND_MAX) - 1.0; // random form -1 to 1
            double phi =
                2.0 * M_PI * ((double)rand() / RAND_MAX); // random 0 to 2pi
            double r = sqrt(1.0 - z * z);
            double x = r * cos(phi);
            double y = r * sin(phi);
            particle_direction.push_back(v3(x, y, z));
        }
    }

    void generate_h_field() {// i think this works
        particle_h_field.resize(params.numberOfParticles, v3(0.0, 0.0, 0.0));
    }

    void generate_e_field() {
        // TODO
    }

    void generate_vel() { 
        particle_velocity.resize(params.numberOfParticles, v3(0.0, 0.0, 0.0));
    }

    void generate_direc_vel() {
        particle_direction_velocity.resize(params.numberOfParticles,
                                           v3(0.0, 0.0, 0.0));
    }

    void first_update_e_dipole() {
        particle_e_dipol.resize(params.numberOfParticles, v3(0.0, 0.0, 0.0));
        update_e_dipole();
    }

    void update_h_field() {
        for (int i = 0; i < params.numberOfParticles; i++) {
            v3 this_h_field = v3(0, 0, 0);
            double prefactor = params.totalMagDipoleMomentParticle / (4 * M_PI);
            for (int j = 0; j < params.numberOfParticles; j++) {
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

    void update_e_field() {
        double prefactor =
            1 / (4 * M_PI * params.relPermittivityMatrix * EPSILON_0);

        for (int i = 0; i < params.numberOfParticles; i++) {
            v3 this_e_field = v3(0, 0, 0);
            for (int j = 0; j < params.numberOfParticles; j++) {
                if (j == i) {
                    continue;
                }
                v3 r_ji = (particle_pos[i] - particle_pos[j]) %
                          params.lengthSimulationCube; // distance with PBC
                v3 r_ji_hat = r_ji.get_direction();
                this_e_field =
                    this_e_field +
                    (prefactor / (pow(r_ji.get_length(), 3))) *
                        (3.0 *
                             (r_ji_hat * (particle_e_dipol[j].dot(r_ji_hat))) -
                         particle_e_dipol[j]);
            }
            particle_e_field[i] = this_e_field;
        }
    }

    void update_e_dipole() {
        double prefactor = params.volumeParticle * EPSILON_0;
        double chi_diff = params.chiEffShortAxisC - params.chiEffLongAxesAB;
        for (int i = 0; i < params.numberOfParticles; i++) {
            v3 left_term = params.chiEffLongAxesAB * particle_e_field[i];
            v3 right_term = chi_diff *
                            (particle_direction[i].dot(particle_e_field[i])) *
                            particle_direction[i];
            particle_e_dipol[i] = prefactor * (left_term + right_term);
        }
    }

    void update_pos_velocity() {
        double h_prefactor =
            3 * MU_0 * pow(params.totalMagDipoleMomentParticle, 2) / (4 * M_PI);
        double e_prefactor =
            3 / (EPSILON_0 * params.relPermittivityMatrix * 2 * M_PI);
        double r_prefactor = 3 * MU_0 *
                             pow(params.totalMagDipoleMomentParticle, 2) /
                             (2 * M_PI * pow(2 * params.eqRadius, 4));
        for (int i = 0; i < params.numberOfParticles; i++) {
            v3 this_h_force = v3(0, 0, 0);
            v3 this_e_force = v3(0, 0, 0);
            v3 this_r_force = v3(0, 0, 0);
            for (int j = 0; j < params.numberOfParticles; j++) {
                if (j == i) {
                    continue;
                }
                v3 r_ji = (particle_pos[i] - particle_pos[j]) %
                          params.lengthSimulationCube; // distance with PBC
                v3 r_ji_hat = r_ji.get_direction();
                this_h_force =
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
                            ((r_ji.get_length() / (2 * params.eqRadius)) - 1)) *
                        r_ji_hat;
            }
            particle_velocity[i] = 1 / params.dragCoeffTransl *
                                   (this_e_force + this_h_force + this_r_force);
        }
    }

    void update_direc_velocity() {
        double h_prefactor = MU_0 * params.totalMagDipoleMomentParticle;
        double e_prefactor =
            params.volumeParticle * EPSILON_0 *
            (params.chiEffShortAxisC - params.chiEffLongAxesAB);
        for (int i = 0; i < params.numberOfParticles; i++) {
            v3 h_cross_d = h_prefactor *
                           (particle_h_field[i] -
                            (particle_direction[i] *
                             (particle_h_field[i].dot(particle_direction[i]))));
            v3 e_cross_d = e_prefactor *
                           (particle_direction[i].dot(particle_e_field[i])) *
                           (particle_e_field[i] -
                            (particle_direction[i] *
                             (particle_e_field[i].dot(particle_direction[i]))));
            particle_direction_velocity[i] =
                1 / params.dragCoeffRot * (h_cross_d + e_cross_d);
        }
    }

    double find_delta_t() {
        double max_value = 0.0;
        for (auto elem : particle_velocity) {
            double norm_velocity = elem.get_length();

            if (norm_velocity > max_value) {
                max_value = norm_velocity;
            }
        }
        double delta_t =
            params.corrFactorVelocity * params.eqRadius / max_value;

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
            particle_direction[i] = (particle_direction[i] +
                                     delta_t * particle_direction_velocity[i])
                                        .get_direction();
        }
    }

    // returns the delta t used for the function
    double take_step() {
        for (int i = 0; i < EL_CONVERGE_ITERATIONS; i++) {
            update_e_field();
            update_e_dipole();
        }
        update_h_field();

        update_direc_velocity();
        update_pos_velocity();

        double delta_t = find_delta_t();

        update_direction(delta_t);
        update_position(delta_t);
        return delta_t;
    }

    void
    run_simulation() { // should we also log the step amount and the
                       // corresponding time bc they are not equally spaced?
        double current_time = 0.0;
        while (current_time < params.simulationTime) {
            current_time += take_step();
        }
    }
};

int main() {
    SimBuilder builder;
    SimulationParameters params = builder.build();

    Sim simulation(params);
    std::cout << "Starting simulation..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    simulation.run_simulation();
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> diff = end - start;

    std::cout << "Simulation finished." << std::endl;
    std::cout << "Time taken: " << diff.count() << " seconds" << std::endl;

    return 0;
}
