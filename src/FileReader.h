/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include "Particle.h"
#include <list>
#include "particleGenerator.h"
#include "utils/Vector.h"
#include "SettingsContainer.h"

class FileReader {

public:
  FileReader();
  virtual ~FileReader();

  void readFile(std::list<Particle> &particles, char *filename);
  void readFile2 (std::list<particleGenerator> &a,std::list<utils::Vector<double,3>> &b,double &end_time, double &delta_t, char* filename);
  SettingsContainer readXMLFile(std::list<particleGenerator> &a, char *filename);
};
