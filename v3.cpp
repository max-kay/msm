#include <cmath>
#include <iostream>
using namespace std;

class v3 {
  public:
    // member variables
    double x, y, z;

    // constructor
    v3(double x, double y, double z) : x(x), y(y), z(z) {}

    // Vector-Scalar Operations
    v3 operator*(const double &num) {
        return v3(x * num, y * num, z * num);
    } // Multiplication
    v3 operator/(const double &num) // Division
    {
        if (num == 0) {
            throw std::domain_error("v3 ERROR: vector-scalar division by zero");
        }
        return v3(x / num, y / num, z / num);
    }

    // Vector-Vector Operations
    v3 operator+(const v3 &vector) {
        return v3(x + vector.x, y + vector.y, z + vector.z);
    } // Addition
    v3 operator-(const v3 &vector) {
        return v3(x - vector.x, y - vector.y, z - vector.z);
    } // Subtraction

    //double operator*(const v3 &vector) {
    //    return double(x * vector.x + y * vector.y + z * vector.z);
    //} // Scalar Product
    double dot(const v3 &vector) {
        return double(x * vector.x + y * vector.y + z * vector.z);
    } // Scalar Product

    // Self-Targeting Operations
    double get_length() {
        return double(std::sqrt(x * x + y * y + z * z));
    } // Euclidian Norm
    v3 get_direction() { return *this / (*this).get_length(); } // Normalize
    v3 operator%(const double &size) // Return to Central Probing Volume
    {
        x -= std::round(x / size) * size;
        y -= std::round(y / size) * size;
        z -= std::round(z / size) * size;

        return v3(x, y, z);
    }
    
};
// Commutative Multiplication
v3 operator*(const double &num, const v3 &vec) {
    return v3(vec.x * num, vec.y * num, vec.z * num);
}