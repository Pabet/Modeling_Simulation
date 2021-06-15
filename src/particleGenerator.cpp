//
// Created by nastya on 07.11.18.
//
#include "particleGenerator.h"
#include <list>
#include "Particle.h"
#include "utils/Vector.h"
#include "thermostat.h"
#include "MaxwellBoltzmannDistribution.h"
#include <vector>


  particleGenerator::particleGenerator(utils::Vector<double,3> c, utils::Vector<double,3> d,float h, float m,utils::Vector<double,3> v,double t,double e,double si)
  {
    ///@brief constructor for the particle Generator in the case that it must generate a cuboid
    ///@param c represents the upper left front corner of the cuboid
    ///@param d represents the number of particles on each dimension, thus the cube having d[0] x d[1] x d[2] particles
        coordinates.operator=(c);
        dimensions.operator=(d);
        mesh_width=h;
        mass=m;
        velocity.operator=(v);
        sigma=si;
        epsilon=e;
        type=t;
    }
 particleGenerator::particleGenerator(utils::Vector<double,3> c, double r,float h, float m,utils::Vector<double,3> v,double t,double e,double si)
 {
    ///@brief constructor for the particle Generator in the case that it must generate a sphere or circle
    ///@param c represents the center of the sphere/circle
    ///@param r represents the radius of the object

        coordinates.operator=(c);
        ///@brief the maximum number of particles in each dimension is 2 * radius / mesh size
        dimensions.operator=(ceil((r*2)/h));
        radius=r;
        isCircle=true;
        mesh_width=h;
        mass=m;
        velocity.operator=(v);
        sigma=si;
        epsilon=e;
        type=t;
    }
void particleGenerator::createGrid(std::list<Particle> &particles, double factor, double target_temp, int dim,bool init_temp,bool membrane, std::vector <Particle> &particles2) {
    ///@brief function for creation of the grid. The  particles in the grid are then saved in the Particle list particles
        for (int i=0;i<dimensions[0];i++){
            double a1=coordinates[0]+mesh_width*i;
            for (int j=0; j<dimensions[1];j++){
                double a2=coordinates[1]+mesh_width*j;

                for (int k=0;k<dimensions[2];k++) {
                    double a [3]={a1,a2,coordinates[2]+mesh_width*k};
                    utils::Vector<double,3>c2 (a);
                    if (!membrane) {
                    Particle p (c2,velocity,mass,type,epsilon,sigma);
                    MaxwellBoltzmannDistribution(p,factor,2);
                    particles.push_back(p);}
                    else {
                        Particle* rightNeighbor=NULL;
                        Particle* rightUpperNeighbor=NULL;
                        Particle* rightLowerNeighbor=NULL;
                        if (i+1<dimensions[0])
                        { rightNeighbor=&particles2[dimensions[1]*(i+1)+j+1];
                        if (j+1<dimensions[1])
                        rightUpperNeighbor=&particles2[dimensions[1]*(i+1)+j+2];
                        if ((j-1)>=0)
                            rightLowerNeighbor=&particles2[dimensions[1]*(i+1)+j];

                        }
                        Particle* UpperBound=NULL;
                        if (j+1<dimensions[1])
                            UpperBound=&particles2[dimensions[1]*i+j+2];
                        Particle* LowerBound=NULL;
                        if (j-1>=0)
                            LowerBound=&particles2[dimensions[1]*i+j];
                        Particle* leftNeighbor=NULL;
                        Particle* leftUpperNeighbor=NULL;
                        Particle* leftLowerNeighbor=NULL;
                        if (i-1>=0){
                            leftNeighbor=&particles2[dimensions[1]*(i-1)+j+1];
                        if (j+1<dimensions[1])
                            leftUpperNeighbor=&particles2[dimensions[1]*(i-1)+j+2];
                        if (j-1>=0)
                            leftLowerNeighbor=&particles2[dimensions[1]*(i-1)+j];
                        }
                           utils::Vector<Particle*,4> neighbors;
                           neighbors[0]=rightNeighbor;
                           neighbors[1]=leftNeighbor;
                           neighbors[2]=UpperBound;
                           neighbors[3]=LowerBound;
                           utils::Vector<Particle*,4> DiagNeighbors;
                           DiagNeighbors[0]=rightUpperNeighbor;
                           DiagNeighbors[1]=rightLowerNeighbor;
                           DiagNeighbors[2]=leftUpperNeighbor;
                           DiagNeighbors[3]=leftLowerNeighbor;
                           Particle p(c2,velocity,mass,neighbors,DiagNeighbors,type, epsilon,sigma);
                        MaxwellBoltzmannDistribution(p,factor,2);
                        particles2.push_back(p);



                    }
                }
        }
    }
   if (init_temp)
    Target_Temperature (target_temp,particles,dim);
}
void particleGenerator::createCircle(std::list<Particle> &particles, double factor ,int dimension,double target_temp,bool init_temp)
{///@brief function for creation of a circle/sphere (based on the @param dimension)
///@brief function first creates a cuboid of size 2r starting at {center[0]-r;center[1]-r;0(or center[2]-r if dimension is 3)}
///@brief the cuboid has the size 2*r in x and y and 2r in z only if dimension is 3, else it has  the size 0 in z
///@param dimension is interpreted as 2 for anything except 3
          for (int i=0;i<dimensions[0];i++){
            double a1=coordinates[0]-radius+mesh_width*i;
            for (int j=0; j<dimensions[1];j++){
                double a2=coordinates[1]-radius+mesh_width*j;
                if (dimension ==3)
                for (int k=0;k<dimensions[2];k++) {
                    double a [3]={a1,a2,coordinates[2]-radius+mesh_width*k};
                    utils::Vector<double,3>c2 (a);
                    utils::Vector<double,3> c3 (a);
                    c3.operator=(c3.operator-(coordinates));
                    if(c3.L2Norm()<=radius)
                    {Particle p (c2,velocity,mass,type,epsilon,sigma);
                    MaxwellBoltzmannDistribution(p,factor,2);
                    particles.push_back(p);
                }}
                else{double a [3]={a1,a2,0};
                    utils::Vector<double,3>c2 (a);utils::Vector<double,3> c3 (a);
                    c3.operator=(c3.operator-(coordinates));
                    ///@brief point is actually generated only if the point is in the radius of the intended circle/sphere
                    if(c3.L2Norm()<=radius)
                    {
                    Particle p (c2,velocity,mass,type,epsilon,sigma);
                    MaxwellBoltzmannDistribution(p,factor,2);
                    particles.push_back(p);}}
        }
    }
    if (init_temp)
    Target_Temperature(target_temp,particles, dimension);
}

   /*   int &particleGenerator::multiGrids(std::list<Particle> &a, double factor) {
        std::list<Particle> b;
                int k=createGrid(b,factor);
                a.splice(a.end(),b);


    }*/


