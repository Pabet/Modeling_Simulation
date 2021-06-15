#pragma once

#include "Particle.h"
#include <list>
class Calculator{
public:
    Calculator();
    virtual ~Calculator();
    void calculateF(std::list<Particle> &particles, bool MBD, double epsilon, double sigma);
    void calculateV(std::list<Particle> &particles,double delta_t);
    void calculateV1P(Particle &p, double delta_t);
    void calculateX1P(Particle &p, double delta_t);
    void calculateX(std::list<Particle> &particles, double delta_t);
};