/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include "utils/Vector.h"

class Particle {

private:
  /** the position of the particle */
  utils::Vector<double, 3> x;

  /** the velocity of the particle */
  utils::Vector<double, 3> v;

  /** the force effective on this particle */
  utils::Vector<double, 3> f;

  /** the force wich was effective on this particle */
  utils::Vector<double, 3> old_f;
utils::Vector<double,3> old_x;
utils::Vector<int,3> PeriodicMovement;
  double epsilon;
  double sigma;
    utils::Vector<Particle*,4> neighbors;
    utils::Vector<Particle*,4> DiagNeighbors;
  /** the mass of this particle */
  double m;

  /** type of the particle. Use it for whatever you want (e.g. to seperate
   * molecules belonging to different bodies, matters, and so on)
   */
  int type;

public:
  Particle(int type = 0);

  Particle(const Particle &other);

  Particle(
      // for visualization, we need always 3 coordinates
      // -> in case of 2d, we use only the first and the second
      utils::Vector<double, 3> x_arg, utils::Vector<double, 3> v_arg,
      double m_arg, int type, double e, double s);
    Particle(
            // for visualization, we need always 3 coordinates
            // -> in case of 2d, we use only the first and the second
            utils::Vector<double, 3> x_arg, utils::Vector<double, 3> v_arg,
            double m_arg,utils::Vector<Particle*,4> n,utils::Vector<Particle*,4> dn, int type , double e, double s);
  virtual ~Particle();

  utils::Vector<double, 3> &getX();

  double &geteps();

  double &getsig();
utils::Vector<Particle*,4> &getneighbors();
utils::Vector<Particle*,4> &getDiagNeighbors();
  utils::Vector<double, 3> &getF();
utils::Vector<double,3> &getOldX();
utils::Vector<int,3> &getPeriodicMovement();
  utils::Vector<double, 3> &getOldF();

  utils::Vector<double, 3> &getV();

  double getM();

  int getType();

  bool operator==(Particle &other);

  std::string toString();
};

std::ostream &operator<<(std::ostream &stream, Particle &p);
