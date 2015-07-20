// Contains templates for 2D and 3D vectors (Vector2<>, Vector3<>) and
// implements common vector operations for them.
// Defines types Vec2d, Vec2i, Vec3d, Vec3d for easy template access.

#pragma once

#include <cmath>

template<typename T>
class Vector2 {
 public:
    Vector2(): x(0), y(0) {};
    Vector2(const T& val): x(val), y(val) {};
    Vector2(const T& xv, const T& yv):
        x(xv), y(yv) {};
    ~Vector2() {};

    double lensqr() const {
        return (pow(x, 2) + pow(y, 2));
    }

    double len() const {
        return sqrt(pow(x, 2) + pow(y, 2));
    }

    bool operator==(const Vector2& other) const {
        return x == other.x && y == other.y;
    }

    bool operator!=(const Vector2& other) const {
        return x != other.x || y != other.y;
    }

    Vector2 operator-(const Vector2& other) const {
        return Vector2(x - other.x, y - other.y);
    }

    Vector2& operator-=(const Vector2& other) {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    Vector2 operator+(const Vector2& other) const {
        return Vector2(x + other.x, y + other.y);
    }

    Vector2& operator+=(const Vector2& other) {
        x += other.x;
        y += other.y;
        return *this;
    }

    template<typename T2>
    Vector2 operator*(const T2& a) const {
        return Vector2(x*a, y*a);
    }

    template<typename T2>
    Vector2& operator*=(const T2& a) {
        x *= a;
        y *= a;
        return *this;
    }

    template<typename T2>
    Vector2 operator/(const T2& a) const {
        return Vector2(x/a, y/a);
    }

    template<typename T2>
    Vector2& operator/=(const T2& a) {
        x /= a;
        y /= a;
        return *this;
    }

    T x, y;
};

template<typename T>
class Vector3 {
 public:
    Vector3(): x(0), y(0), z(0) {};
    Vector3(const T& val): x(val), y(val), z(val) {};
    Vector3(const T& xv, const T& yv, const T& zv):
        x(xv), y(yv), z(zv) {};
    Vector3(const Vector2<T>& v2, const T& zv):
        x(v2.x), y(v2.y), z(zv) {};
    ~Vector3() {};

    Vector3& zero() {
        x = 0;
        y = 0;
        z = 0;
        return *this;
    }

    T lensqr() const {
        return (pow(x, 2) + pow(y, 2) + pow(z, 2));
    }

    double len() const {
        return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    }

    Vector3 normalized() const {
        double l = len();
        if (l > 1e-10) {
            return Vector3(x / l, y / l, z / l);
        } else {
            return Vector3(0);
        }
    }

    T dot(const Vector3& other) const {
        return (x * other.x + y * other.y +  z * other.z);
    }

    Vector3 cross(const Vector3& other) const {
        return Vector3(y*other.z - z*other.y, z*other.x - x*other.z, x*other.y - y*other.x);
    }

    Vector2<T> getXY() const {
        return Vector2<T>(x, y);
    }

    bool operator==(const Vector3& other) const {
        return x == other.x && y == other.y && z == other.z;
    }

    bool operator!=(const Vector3& other) const {
        return x != other.x || y != other.y || z != other.z;
    }

    Vector3 operator-(const Vector3& other) const {
        return Vector3(x - other.x, y - other.y, z - other.z);
    }

    Vector3& operator-=(const Vector3& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    Vector3 operator+(const Vector3& other) const {
        return Vector3(x + other.x, y + other.y, z + other.z);
    }

    Vector3& operator+=(const Vector3& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    template<typename T2>
    Vector3 operator*(const T2& a) const {
        return Vector3(x*a, y*a, z*a);
    }

    template<typename T2>
    Vector3& operator*=(const T2& a) {
        x *= a;
        y *= a;
        z *= a;
        return *this;
    }

    template<typename T2>
    Vector3 operator/(const T2& a) const {
        return Vector3(x/a, y/a, z/a);
    }

    template<typename T2>
    Vector3& operator/=(const T2& a) {
        x /= a;
        y /= a;
        z /= a;
        return *this;
    }

    T x, y, z;
};

typedef Vector2<double> Vec2d ;
typedef Vector2<int> Vec2i;
typedef Vector3<double> Vec3d ;
typedef Vector3<int> Vec3i;

template<typename T1, typename T2>
Vector3<T2> operator*(const T1& a, const Vector3<T2> vec){
    return Vector3<T2>(a * vec.x, a * vec.y, a * vec.z);
}

template<typename T1, typename T2>
Vector2<T2> operator*(const T1& a, const Vector2<T2> vec){
    return Vector2<T2>(a * vec.x, a * vec.y);
}
