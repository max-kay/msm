#include "builder.cpp"
#include <iostream>
#include <vector>
#define EPSILON_0 (double) 8.8541878188e-12
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

    void update_e_field() {
        double prefactor =
            1 /
            (4 * M_PI *
             params.relPermittivityParticle*EPSILON_0); // NO DIELECTRIC CONST. epslion_0

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
                        (3.0 *
                             (r_ji_hat * (particle_e_dipol[j].dot(r_ji_hat))) -
                         particle_e_dipol[j]);
            }
            particle_e_field[i] = this_e_field; // update the list
        }
    }
    void update_h_field() {
        // TODO
    }
};

int main() {
    std::cout << "Hello World" << std::endl;
    return 0;
}
