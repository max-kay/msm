#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <vector>
#include <string>
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

int write_numpy_file(const std::vector<v3> &buf, const std::string &path) {
    size_t len = buf.size();
    FILE *f = fopen(path.c_str(), "wb");
    if (!f)
        return -1;

    /* Build header dict string */
    char dict[256];
    int dict_len = snprintf(
        dict, sizeof(dict),
        "{'descr': '<f8', 'fortran_order': False, 'shape': (%zu, 3), }", len);
    if (dict_len < 0 || dict_len >= (int)sizeof(dict)) {
        fclose(f);
        return -1;
    }

    /*
     * Total = 10 (magic + version + header_len u16) + HEADER_LEN
     * HEADER_LEN must make the total a multiple of 64.
     * The header ends with '\n'; everything before is space-padded.
     */
    const int PREAMBLE = 12;
    int total_unpadded = PREAMBLE + dict_len + 1; /* +1 for '\n' */
    int remainder = total_unpadded % 64;
    int padding = (remainder == 0) ? 0 : (64 - remainder);
    int header_len = dict_len + padding + 1; /* dict + spaces + newline */

    if (header_len > 65535) {
        fclose(f);
        return -1;
    }

    uint32_t hlen_le = (uint32_t)header_len;
    /* Write magic + version */
    if (fwrite("\x93NUMPY\x02\x00", 1, 8, f) != 8)
        goto err;

#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    hlen_le = __builtin_bswap32(hlen_le);
#endif
    if (fwrite(&hlen_le, 4, 1, f) != 1) /* 4 bytes, not 2 */
        goto err;

    /* Write dict, space padding, newline */
    if (fwrite(dict, 1, (size_t)dict_len, f) != (size_t)dict_len)
        goto err;
    for (int i = 0; i < padding; i++)
        if (fputc(' ', f) == EOF)
            goto err;
    if (fputc('\n', f) == EOF)
        goto err;

    for (size_t i = 0; i < len; i++) {
        double xyz[3] = {buf[i].x, buf[i].y, buf[i].z};
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
        for (int j = 0; j < 3; j++) {
            uint64_t tmp;
            memcpy(&tmp, &xyz[j], 8);
            tmp = __builtin_bswap64(tmp);
            memcpy(&xyz[j], &tmp, 8);
        }
#endif
        if (fwrite(xyz, 8, 3, f) != 3)
            goto err;
    }

    fclose(f);
    return 0;

err:
    fclose(f);
    return -1;
}
