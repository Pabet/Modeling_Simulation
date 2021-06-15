#pragma once

#include <list>
#include "utils/Vector.h"
#include "particleContainer.h"
#include "Particle.h"
#include "Calculator.h"
class Mesh{
public:
    std::list<Particle> particles;
    double gravitationa;
    particleContainer* Cells;
    double delta_t;
    double t_target;
    bool inittemp;
    int xsize;
    int ysize;
    int zsize;
    bool membrane=false;
    double k=0.0;
    double r0=0.0;
    double F_z_UP=0.0;
    double rl=0.0;
    double cutoff;
    utils::Vector<double, 3> domain_size;
    bool zdim;
    int x1boundary=0; ///@brief boundary rules: 0= no boundary, especially useful in the case of 2 dimensions for the third dimension
    int x2boundary=0;/// 1=outflow ;2=reflecting ;3 = periodic;
    int y1boundary=0;
    int y2boundary=0;
    int z1boundary=0;
    int z2boundary=0;
    int parallelisationMethod=0;
    Calculator s;
    void Target_Temperature (double t_target);
    Mesh(bool inittemp,double t_target,double rcutoff,utils::Vector<double, 3> domain_size,std::list <Particle> &tempparticles,double d_t,double gravtiation,int x1boundary,int x2boundary,int y1boundary,int y2boundary,int z1boundary,int z2boundary, int parallelisationMethod,double r2);
     Mesh(bool inittemp,double t_target, double rcutoff,utils::Vector<double, 3> domain_size,std::vector<Particle>&tempparticles,double d_t,double gravitation,int x1bound,int x2bound,int y1bound,int y2bound,int z1bound,int z2bound,double k1,double r1,double F_up,double r2,int pMethod);
    void calculateV();
   double getRDF();
    void reassignMesh();
    void calculateX(bool rdf);
    void calculateF(bool timem);
    void halloCellsDecisions();
};
