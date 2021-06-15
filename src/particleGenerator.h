//
// Created by nastya on 08.11.18.
//

#pragma once
#include <list>
#include "Particle.h"
#include "utils/Vector.h"
#include "MaxwellBoltzmannDistribution.h"
#include <vector>
class particleGenerator{
public:
  double velocityBM;
  utils::Vector<double,3> coordinates;
  utils::Vector<double,3> dimensions;
  float mesh_width;
  float mass;
  bool isCircle=false;
  double radius;
  double epsilon;
  double sigma;
  double type;
  utils::Vector<double,3> velocity;
public:
    particleGenerator(utils::Vector<double,3> c, utils::Vector<double,3> d,float h, float m,utils::Vector<double,3> v,double t,
                      double e, double si);
    particleGenerator(utils::Vector<double,3> c, double r,float h, float m,utils::Vector<double,3> v,double t, double e,
                      double si);


void createGrid (std::list<Particle> &particles, double factor,double target_temp, int dimension, bool init_temp,bool membrane,std::vector <Particle> &particles2);

void createCircle (std::list<Particle> &particles, double factor, int dimension, double target_temp, bool init_temp);
  //  int &multiGrids (std::list<Particle> &a,double factor);

};

