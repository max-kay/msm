#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int main() {
    srand(time(nullptr)); // seed the generator

    double rng_1 = (double) rand() / RAND_MAX-0.5;
    double rng_2 = (double) rand() / RAND_MAX-0.5;
    double rng_3 = (double) rand() / RAND_MAX-0.5;

    cout << rng_1 << '\n';
    cout << rng_2 << '\n';
    cout << rng_3 << '\n';

    return 0;
}