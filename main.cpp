#include "SimulationParameters_and_SimBuilder.cpp"
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

// corrFactorVelocity
class Sim {
  public:
    SimulationParameters params;
    vector<v3>
        particle_e_field; // platz wird reserviert mit einem typ -- einzonen
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

    void update_e_field() { // funktion ohne stützräder (max) geschrieben -->
                            // überprüfen SEHR nötig
        double prefactor =
            1 /
            (4 * M_PI *
             params.relPermittivityParticle); // NO DIELECTRIC CONST. epslion_0
                                              // bc we are not using it right?
        for (int i = 0; params.numberOfParticles; i++) {
            vector<v3> dist_to_particle_i;
            for (int j = 0; params.numberOfParticles; i++) {
                dist_to_particle_i[j] =
                    particle_pos[i] -
                    particle_pos[j]; // richtig herum? fine with periodic bc?
                // these are here to calculate the distance from our i-th
                // particle to the j-th particle

                vector<v3> e_field_change_of_i;
                for (int k = 0; params.numberOfParticles; i++) {
                    if (k != i) { // dont devide by 0
                        e_field_change_of_i[k] =
                            prefactor /
                            (pow(dist_to_particle_i[k].get_length(),
                                 3)) *
                            3 *
                            ((particle_e_dipol[j].dot(
                                  dist_to_particle_i[k].get_direction()))
                                 .dot(dist_to_particle_i[k].get_direction()) -
                             particle_e_dipol[k]); // wrong implematation of dot
                                                   // functio idk how to use
                                                   // correctly
                    }
                }
            }
        }
    }
};

int main() {
    std::cout << "Hello World" << std::endl;
    return 0;
}
