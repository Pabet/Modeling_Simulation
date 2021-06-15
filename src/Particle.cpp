/*
* Particle.cpp
        *
        *  Created on: 23.02.2010
*      Author: eckhardw
*/

#include "Particle.h"
#include <iostream>
#include <sstream>

/* Iterator Pattern: hilft beim wenn man verschiedene Objekte iterieren möchte
 * bei Atomen/Molekülen wäre auch das Composite Pattern denkbar */
template <class Item>
class Iterator {
public:
    virtual void First() = 0;
    virtual void Next() = 0;
    virtual bool IsDone() const = 0;
    virtual Item CurrentItem() const = 0;
protected:
    Iterator();
};


class ParticleIterator : public Iterator<Particle> {
public:
    /*hier könnte man dann einen speziellen Iterator
     *angeben, der Moleküle oder Partikel durchgeht */
};

Particle::Particle(int type_arg) {
    type = type_arg;
  //  std::cout << "Particle generated!" << std::endl;
    f = 0.0;
    old_f = 0.0;

}

Particle::Particle(const Particle &other) {
    x = other.x;
    v = other.v;
    f = other.f;
    old_f = other.old_f;
    m = other.m;
    type = other.type;
    epsilon=other.epsilon;
    sigma=other.sigma;
    neighbors.operator=(other.neighbors);
    DiagNeighbors.operator=(other.DiagNeighbors);
    PeriodicMovement[0]=0;
    PeriodicMovement[1]=0;
    PeriodicMovement[2]=0;
    old_x=other.old_x;
   // std::cout << "Particle generated by copy!" << std::endl;
    ;
}

// Todo: maybe use initializater list instead of copy?
Particle::Particle(utils::Vector<double, 3> x_arg,
                   utils::Vector<double, 3> v_arg, double m_arg, int type_arg,double e, double s) {
    x = x_arg;
    v = v_arg;
    m = m_arg;
    type = type_arg;
    f = 0.0;
    old_f = 0.0;
    epsilon=e;
    sigma=s;
    neighbors[0]=NULL;
    neighbors[1]=NULL;
    neighbors[2]=NULL;
    neighbors[3]=NULL;
    DiagNeighbors[0]=NULL;
    DiagNeighbors[1]=NULL;
    DiagNeighbors[2]=NULL;
    DiagNeighbors[3]=NULL;
    PeriodicMovement[0]=0;
    PeriodicMovement[1]=0;
    PeriodicMovement[2]=0;
    old_x.operator=(x);
   /// std::cout << "Particle generated!" << std::endl;
}
Particle::Particle(utils::Vector<double, 3> x_arg,
                   utils::Vector<double, 3> v_arg, double m_arg,utils::Vector<Particle*,4> n,utils::Vector<Particle*,4> dn, int type_arg,double e, double s) {
    x = x_arg;
    v = v_arg;
    m = m_arg;
    type = type_arg;
    f = 0.0;
    old_f = 0.0;
    epsilon=e;
    sigma=s;
    neighbors.operator=(n);
    DiagNeighbors.operator=(dn);
    PeriodicMovement[0]=0;
    PeriodicMovement[1]=0;
    PeriodicMovement[2]=0;
    old_x.operator=(x);
    /// std::cout << "Particle generated!" << std::endl;
}
Particle::~Particle() {
    /// std::cout << "Particle destructed!" << std::endl;
    }

utils::Vector<double, 3> &Particle::getX() { return x; }
double &Particle::geteps() { return epsilon;}
double &Particle::getsig() {return sigma;}

utils::Vector<double, 3> &Particle::getV() { return v; }

utils::Vector<double, 3> &Particle::getF() { return f; }
utils::Vector<double ,3> &Particle::getOldX() {return old_x;}
utils::Vector<int,3> &Particle::getPeriodicMovement() {return PeriodicMovement;}
utils::Vector<double, 3> &Particle::getOldF() { return old_f; }
utils::Vector<Particle*,4> &Particle::getneighbors() { return neighbors;}
utils::Vector<Particle*,4> &Particle::getDiagNeighbors() {return DiagNeighbors;}
double Particle::getM() { return m; }

int Particle::getType() { return type; }

std::string Particle::toString() {
    std::stringstream stream;
    stream << "Particle: X:" << x << " v: " << v << " f: " << f
           << " old_f: " << old_f <<"m: "<<m <<" type: " << type<<" epsilon: "<<epsilon<<" sigma: "<<sigma;
    return stream.str();
}

bool Particle::operator==(Particle &other) {
    if ((x == other.x) && (v == other.v) && (f == other.f) &&
        (type == other.type) & (m == other.m) && (old_f == other.old_f)&&(sigma==other.sigma)&&(epsilon==other.epsilon)) {
        return true;
    }

    return false;
}

std::ostream &operator<<(std::ostream &stream, Particle &p) {
    stream << p.toString();
    return stream;
}