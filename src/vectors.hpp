#pragma once

template<typename T>
class Vector3 {
 public:
    Vector3(): x(0), y(0), z(0) {};
    Vector3(T val): x(val), y(val), z(val) {};
    Vector3(T xv, T yv, T zv):
        x(xv), y(yv), z(zv) {};
    ~Vector3() {};
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
    T x, y;
};

class Vec2d: public Vector2<double> {
};

class Vec2i: public Vector2<int> {
};
