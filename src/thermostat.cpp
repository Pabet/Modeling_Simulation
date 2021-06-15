//
// Created by nastya on 11.12.18.
//

#include <list>
#include "thermostat.h"
#include "Particle.h"
#include "utils/Vector.h"
using namespace std;


void Target_Temperature (double t_target,  std::list<Particle> &particles,int dimensions){


    std::list<Particle>::iterator iterator=particles.begin();
    double e_kin=0;
    while (iterator!=particles.end()){
        Particle p=*iterator;
       if (dimensions==2)
        e_kin=e_kin+0.5*(p.getV()[0]*p.getV()[0]+p.getV()[1]*p.getV()[1]);
        if (dimensions==3)
            e_kin=e_kin+0.5*(p.getV()[0]*p.getV()[0]+p.getV()[1]*p.getV()[1]+p.getV()[2]*p.getV()[2]);
        iterator++;}
    double temp=2*e_kin/(dimensions*particles.size());

    double betha=sqrt(t_target/temp);

 //  if (betha!=1) {
     std::list<Particle>::iterator  iterator1 = particles.begin();
       while (iterator1 != particles.end()) {
           Particle &p = *iterator1;


           p.getV().operator=(p.getV().operator*(betha));

           iterator1++;
       }
 //  }




    }

