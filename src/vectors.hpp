#pragma once

#include <cmath>

template<typename T>
class Vector3 {
 public:
    Vector3(): x(0), y(0), z(0) {};
    Vector3(T val): x(val), y(val), z(val) {};
    Vector3(T xv, T yv, T zv):
        x(xv), y(yv), z(zv) {};
    ~Vector3() {};

    T lensqr() const {
        return (pow(x, 2) + pow(y, 2) + pow(z, 2));
    }

    double len() const {
        return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    }

    Vector3 normalized() const {
        double l = len();
        return Vector3(x / l, y / l, z / l);
    }

    T dot(const Vector3& other) const {
        return (x * other.x + y * other.y +  z * other.z);
    }

    Vector3 operator-(const Vector3& other) const {
        return Vector3(x - other.x, y - other.y, z - other.z);
    }

    Vector3 operator+(const Vector3& other) const {
        return Vector3(x + other.x, y + other.y, z + other.z);
    }

    T x, y, z;
};

typedef Vector3<double> Vec3d ;
typedef Vector3<int> Vec3i;

template<typename T>
class Vector2 {
 public:
    Vector2(): x(0), y(0) {};
    Vector2(T& val): x(val), y(val) {};
    Vector2(T& xv, T& yv):
        x(xv), y(yv) {};
    ~Vector2() {};

    Vector2 operator-(const Vector2& other) const {
        return Vector2(x - other.x, y - other.y);
    }

    Vector2 operator+(const Vector2& other) const {
        return Vector2(x + other.x, y + other.y);
    }

    double lensqr() const {
        return (pow(x, 2) + pow(y, 2));
    }

    double len() const {
        return sqrt(pow(x, 2) + pow(y, 2));
    }

    T x, y;
};

class Vec2d: public Vector2<double> {
};

class Vec2i: public Vector2<int> {
};
