//
// Created by nastya on 30.11.18.
//
#pragma once
#include "utils/Vector.h"
#include <list>
#include "Particle.h"
#include <vector>
class particleGenerator2{
public:
    double radius;
    bool voutflow;
    bool vperiodic;
    bool vreflecting;
    bool houtflow;
    bool hperiodic;
    bool hreflecting;
    double x;
    double y;
    double z;
    bool mbd;
  

public:
    particleGenerator2(bool vo,bool vp, bool vrf,bool ho,bool hp, bool hrf, utils::Vector<double,3> d, double r,bool f);
    int &createMesh (std::list<Particle> &particles,double sigma, double epsilon,bool gravity,double g);

};

