//
// Created by Victor Stroescu on 14/11/2018.
//
#include "Calculator.h"
#include <list>
#include "outputWriter/VTKWriter.h"
//#include "MaxwellBoltzmannDistribution.h"
#include "particleGenerator.h"

using namespace std;

Calculator::Calculator() {}

Calculator::~Calculator() {}

void Calculator::calculateF(std::list<Particle> &particles, bool MBD,  double epsilon,double sigma) {

    ///@brief calculate the Force afecting the particle based on the distance to the other particles and their mass
    list<Particle>::iterator iterator;
    iterator = particles.begin();
    int counter=0;
    double l [3]={0,0,0};
    vector<utils::Vector<double,3>>forces;
    ///@param MBD is true if the Program uses the MaxwellBoltzmanDistribution, else false
    ///@brief this is done in order to keep the functionality of task 1
    if(MBD) {
        for (long unsigned int i = 0; i < particles.size() - 1; ++i) {
            utils::Vector<double, 3> a(l);
            forces.push_back(a);
        }
    }
    ///@brief for each particle, determined by the iterator, calculate all the forces that affect the particle aand summ them up
    while (iterator != particles.end()) {

        list<Particle>::iterator innerIterator;
        if(MBD) {
            innerIterator = iterator;

            innerIterator++;

            utils::Vector<double, 3> c(l);
            Particle &p1 = *iterator;

            int counter1 = counter;
            while (innerIterator != particles.end()) {

                Particle &p2 = *innerIterator;
                float k = (p1.getX() - p2.getX()).L2Norm();
                float l3 = sqrt(p1.getsig()*p2.getsig()) / k;
                float f1 = ((24 * ((p1.geteps()+p2.geteps())/2)) / pow(k, 2)) * (pow(l3, 6) - 2 * pow(l3, 12));
                utils::Vector<double, 3> b = (p2.getX().operator-(p1.getX())).operator*(f1);
                c.operator=(c.operator+(b));

                forces[counter1].operator=(forces[counter1].operator+(b.operator*(-1)));

                counter1++;
                innerIterator++;
            }
            p1.getOldF().operator=(p1.getF());
            p1.getF().operator=(c);
            counter++;

    }else{
            innerIterator = particles.begin();
            vector<utils::Vector<double, 3>> forces;
            ///@brief Forces affecting Particle i
            while (innerIterator != particles.end()) {
                if (innerIterator != iterator) {
                    Particle &p1 = *iterator;
                    Particle &p2 = *innerIterator;
                    // @TODO: insert calculation of force here!
                    //i:p1, j:p2
                    float scalar = (p1.getM() * p2.getM()) /
                                   pow((p1.getX().operator-(p2.getX())).L2Norm(), 3);      ///@brief scalar
                    utils::Vector<double, 3> force = (p2.getX().operator-(p1.getX())).operator*(scalar);                         ///F(i,j) = vector*scalar
                    forces.push_back(force);                                                                ///@brief pushback F(i,j)
                }
                ++innerIterator;
            }
            utils::Vector<double, 3> newForce;
            newForce.operator=(0.0);
            for (unsigned int i = 0; i < forces.size(); i++) {
                newForce = newForce.operator+(forces[i]);
            }

            Particle &p1 = *iterator;
            p1.getOldF().operator=(p1.getF());
            p1.getF().operator=(newForce);
            ++iterator;

    }
    }

    if(MBD){
        list<Particle>::iterator iterator3;
        iterator3 = particles.begin();

        iterator3++;

        for (long unsigned int i = 0; i < particles.size() - 1; ++i) {
            Particle &p1 = *iterator3;
            p1.getF().operator=((p1.getF().operator+(forces[i])));
            iterator3++;

        }
    }

}
void Calculator::calculateV1P(Particle &p, double delta_t) {
    ///@brief Calculate the velocity based on the Forces

        p.getV().operator= (p.getV().operator+((p.getOldF().operator+(p.getF())).operator*(delta_t/(2*p.getM()))));

}
void Calculator::calculateV(std::list<Particle> &particles, double delta_t) {
    ///@brief Calculate the velocity based on the Forces
    list<Particle>::iterator iterator = particles.begin();
    while (iterator != particles.end()) {
        Particle &p = *iterator;
        p.getV().operator= (p.getV().operator+((p.getOldF().operator+(p.getF())).operator*(delta_t/(2*p.getM()))));
        ++iterator;

    }
}
void Calculator::calculateX1P(Particle &p, double delta_t) {
    double x_new [] ={0,0,0};
    for(int i=0;i<3;i++)
    {
        x_new[i]=p.getX().operator[](i)+delta_t*p.getV().operator[](i)+delta_t*delta_t*p.getF().operator[](i)/(2*p.getM());
    }
    p.getX().operator=(utils::Vector<double,3> (x_new) );
}
void Calculator::calculateX(std::list<Particle> &particles, double delta_t) {
    ///@brief calculate the new position
    list<Particle>::iterator iterator = particles.begin();
    while (iterator != particles.end()) {

        Particle &p = *iterator;
        double x_new [] ={0,0,0};
        /// @brief For each dimmension calculate the movement based on the velocity of the dimension
        for(int i=0;i<3;i++)
        {
            x_new[i]=p.getX().operator[](i)+delta_t*p.getV().operator[](i)+delta_t*delta_t*p.getF().operator[](i)/(2*p.getM());
        }
        p.getX().operator=(utils::Vector<double,3> (x_new) );
        ++iterator;
    }
}


