//
// Created by Victor Stroescu on 04/01/2019.
//
#pragma once
#include "Particle.h"
#include <vector>
class particleContainer{
private:
    int elementNumber;
    std::vector<Particle*> particleArray;
public:
    particleContainer();
    int getSize();
    void addParticle (Particle* particle);
    void removeparticleAt(int s);
    Particle* getParticleAt(int i);
    int getElementNumber();
};