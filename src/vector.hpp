template<typename T>
class Vector3 {
 public:
    Vector3(): x(0), y(0), z(0) {};
    ~Vector3() {};
    T x, y, z;
};

class Vec3d: public Vector3<double> {
};

class Vec3i: public Vector3<int> {
};

template<typename T>
class Vector2 {
 public:
    Vector2(): x(0), y(0) {};
    ~Vector2() {};
    T x, y;
};

class Vec2d: public Vector2<double> {
};

class Vec2i: public Vector2<int> {
};
