//
// Created by nastya on 18.12.18.
//
//
// Created by nastya on 17.12.18.
//

#pragma once

#include "Particle.h"
#include <list>
void writeCheckpoint (std::list<Particle> particles,const std::string &filename, double iteration);
void readCheckpoint (std::list<Particle> &particles,const std::string filename,  double &iteration);