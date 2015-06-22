#pragma once

#include <memory>

#include <integrators.hpp>

enum MinimiserType {
    STEEPEST_DESCENT,
    FIRE
};

class Simulation;
class System;

class Minimiser {
 public:
    virtual ~Minimiser() {};
    virtual void initialize(Simulation& simulation) = 0;
    virtual void minimise(System& system) = 0;
 private:
    IntegratorType integrator_type_;
    double dt;
};

class SDMinimiser: public Minimiser {
 public:
    SDMinimiser() {};
    ~SDMinimiser() {};
    void initialize(Simulation& simulation) override;
    void minimise(System& system) override;
};

class FIREMinimiser: public Minimiser {
 public:
    FIREMinimiser() {};
    ~FIREMinimiser() {};
    void initialize(Simulation& simulation) override;
    void minimise(System& system) override;
};
