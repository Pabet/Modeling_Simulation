#pragma once

#include "Particle.h"

#include <list>

class Input {
public:
    double end_timei;
    double delta_ti;
    Input() ;
std::list<Particle> InterpretInput(int c,char*param[]);};
