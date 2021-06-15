//
// Created by Victor Stroescu on 04/01/2019.
//

#include "utils/Vector.h"
#include <list>
#include "particleContainer.h"
#include "Particle.h"
#include <memory.h>
#include <stdexcept>
particleContainer::particleContainer() {
    elementNumber=0;
    particleArray=std::vector<Particle*>();
}
int particleContainer::getElementNumber(){ return elementNumber;}
void particleContainer::removeparticleAt(int s)
{
    if(s>=elementNumber)
    { std::cout<<s<<elementNumber;
        throw std::runtime_error("Invalid acces of particle Container"); }
    particleArray.erase(particleArray.begin()+s);
    elementNumber--;
}
void particleContainer::addParticle (Particle* particle)
{
    elementNumber++;;
     particleArray.push_back(particle);
}

Particle* particleContainer::getParticleAt(int i)
{
    if(i<elementNumber)
    {return particleArray[i];}
    else
        throw std::runtime_error("Invalid acces of particle Container");
}
