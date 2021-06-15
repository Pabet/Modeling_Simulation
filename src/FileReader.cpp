/*
 * FileReader.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "FileReader.h"
#include "utils/Vector.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#ifdef logOn
#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"
#include "log4cxx/helpers/exception.h"
using namespace log4cxx;
using namespace log4cxx::helpers;
LoggerPtr fileLogger(Logger::getLogger("MolSim.File"));
#else
#define LOG4CXX_TRACE(logger, expression)
#define LOG4CXX_DEBUG(logger, expression)
#define LOG4CXX_INFO(logger, expression)
#define LOG4CXX_WARN(logger, expression)
#define LOG4CXX_ERROR(logger, expression)
#define LOG4CXX_FATAL(logger, expression)
#endif

#include "shapes.hxx"

using namespace std;

FileReader::FileReader() {}

FileReader::~FileReader() {}


void FileReader::readFile(std::list<Particle> &particles, char *filename) {

   /// @brief FIle reader for first task
  double x[] = {0, 0, 0};
  double v[] = {1, 1, 1};
  double m = 1;
  int num_particles = 0;

  std::ifstream input_file(filename);
  string tmp_string;

  if (input_file.is_open()) {

    getline(input_file, tmp_string);
    LOG4CXX_INFO(fileLogger, "Read line: " << tmp_string);

    while (tmp_string.size() == 0 || tmp_string[0] == '#') {
      getline(input_file, tmp_string);
      LOG4CXX_INFO(fileLogger, "Read line: " << tmp_string);
    }

    istringstream numstream(tmp_string);
    numstream >> num_particles;
    LOG4CXX_INFO(fileLogger, "Reading " << num_particles << ".");
    getline(input_file, tmp_string);
    LOG4CXX_INFO(fileLogger, "Read line: " << tmp_string);

    for (int i = 0; i < num_particles; i++) {
      istringstream datastream(tmp_string);

      for (int j = 0; j < 3; j++) {
        datastream >> x[j];
      }
      for (int j = 0; j < 3; j++) {
        datastream >> v[j];
      }
      if (datastream.eof()) {
        LOG4CXX_ERROR(fileLogger,
            "Error reading file: eof reached unexpectedly reading from line "
            << i);
        exit(-1);
      }
      datastream >> m;
      Particle p(x, v, m);
      particles.push_back(p);

      getline(input_file, tmp_string);
      LOG4CXX_INFO(fileLogger, "Read line: " << tmp_string);
    }
  } else {
    LOG4CXX_ERROR(fileLogger, "Error: could not open file " << filename);
    exit(-1);
  }
}
void FileReader::readFile2 (std::list<particleGenerator> &a, std::list<utils::Vector<double,3>> &b,double &end_time,double &delta_t , char *filename) {
    /// @brief FIle reader for second task
    int s = 0;
    double r = 1;
    double x[] = {0, 0, 0};
    double v[] = {1, 1, 1};
    double N[] = {0, 0, 0};
    double m = 1;                //neu?
    float h = 0;                    //neu
    double epsilon_sigma_factor[3];

    //utils::Vector<double,3> vBM(x);

    int num_cuboids = 0;
    std::ifstream input_file(filename);
    string tmp_string;

    if (input_file.is_open()) {

        getline(input_file, tmp_string);
        LOG4CXX_INFO(fileLogger, "Read line: " << tmp_string);


        while (tmp_string.size() == 0 || tmp_string[0] == '#') {
            getline(input_file, tmp_string);
            LOG4CXX_INFO(fileLogger, "Read line: " << tmp_string);

        }

        istringstream numstream(tmp_string);
        ///@brief num_cuboids the number of objects
        numstream >> num_cuboids;
        LOG4CXX_INFO(fileLogger, "Reading " << num_cuboids << ".");

        getline(input_file, tmp_string);
        LOG4CXX_INFO(fileLogger, "Read line: " << tmp_string);

        /// @brief each of the objects the data is saved depending on the type of object. Each object is on a line of the input file. Number of objects is equal to num_cuboids
        for (int i = 0; i < num_cuboids; i++) {
            istringstream datastream(tmp_string);
            /// @brief if s is 1 then object is a circle else if it is 0 its a cuboid
            datastream>>s;
            ///@brief if circle then x is center else its upper right corner of the cuboid
            for (int j = 0; j < 3; j++) {
                datastream >> x[j];
            }
            ///@brief v is the initial velocity of the object
            for (int j = 0; j < 3; j++) {
                datastream >> v[j];
            }
            ///@brief if the object is a cuboid then the next three values are the number of points in each dimension else it is the radius of the circle
            if (s==0)
                {for (int j=0; j<3;j++){
                datastream>>N[j];}}
                else {
                datastream >> r;
                }

            if (datastream.eof()) {
                LOG4CXX_ERROR(fileLogger,
                              "Error reading file: eof reached unexpectedly reading from line "
                                      << i);

                exit(-1);
            }
            ///@brief m is the mass of the particles of the object
            datastream>>m;
            ///@brief h is the mesh width
            datastream>>h;
            ///@brief epsilon sigma and factor are saved in the array epsilon_sigma_factor at position 0 1 and respectively 2
            datastream>>epsilon_sigma_factor[0];
            datastream>>epsilon_sigma_factor[1];
            datastream>>epsilon_sigma_factor[2];


            double vBM=epsilon_sigma_factor[2]*sqrt(8/3.14); 			//Neu
            ///@brief if the object is a cuboid the constructor for the particleGenerator is called for a cuboid with all the required parameters
            if (s==0){
            particleGenerator p (x,N,h,m,v,vBM,epsilon_sigma_factor[0],epsilon_sigma_factor[1]);
            a.push_back(p);}
            else{
                ///@brief if the object is a sphere the constructor for the particleGenerator is called for a sphere with all the required parameters
                particleGenerator p(x,r,h,m,v,vBM,epsilon_sigma_factor[0],epsilon_sigma_factor[1]);
                a.push_back(p);
            }
            utils::Vector<double, 3> a1(epsilon_sigma_factor);

            b.push_back(a1);
            getline(input_file, tmp_string);
            LOG4CXX_INFO(fileLogger, "Read line: " << tmp_string);}
        ///@brief end of objects

        istringstream datastream1(tmp_string);
        ///@brief last line has end_time and delta_time
        datastream1>>end_time;
        datastream1>>delta_t;
    } else {
        LOG4CXX_ERROR(fileLogger, "Error: could not open file " << filename);

        exit(-1);
    }
}

SettingsContainer FileReader::readXMLFile(std::list<particleGenerator> &a, char *filename){

    double delta_t;
    double end_time;
    double factor;
    utils::Vector<double, 3> domain_size;
    double rcutoff;
    int x1_boundary_condition = 0;
    int x2_boundary_condition = 0;
    int y1_boundary_condition = 0;
    int y2_boundary_condition = 0;
    int z1_boundary_condition = 0;
    int z2_boundary_condition = 0;
    bool brownian_motion=false;
    double initial_temperature;
    double nthermostat;
    double target_temperature;
    double temperature_difference;
    double gravitation;
    double r0 = 0.0;
    double k = 0.0;
    double fz_up = 0.0;
    double rl = 0.0;

    bool writeCheckpoint;
    bool readCheckpoint;
    double writeCheckpointTime;

    int parallelisation_method = 0;
    int force_calculation_method = 0; //0:Lennard-Jones, 1:Smoothed-L-J, 2:gravitation

    /// @brief FIle reader for third task
    double t=0;
    double r=1;
    double e=0;
    double si=0;
    double x[] = {0, 0, 0};
    double v[] = {1, 1, 1};
    double N[]={0,0,0};
    double m;				//neu?
    float  h=0;					//neu

    try {
        ///@brief try to find delta_t, end_time, epsilon, sigma and factor, if they are not found use the standard value
        unique_ptr<shapesType> s1(shapes(filename));

        //settings
        string dt = s1->settings().delta_t().t().get();
        delta_t = std::stod(dt);
        LOG4CXX_INFO(fileLogger, "delta_t: " << dt);

        string et = s1->settings().end_time().t().get();
        end_time = std::stod(et);
        LOG4CXX_INFO(fileLogger, "end_time: " << et);

        string ft = s1->settings().factor().val().get();
        factor = std::stod(ft);
        LOG4CXX_INFO(fileLogger, "factor: " << ft);

        string dx = s1->settings().domain_size().x().get();
        domain_size[0] = std::stod(dx);
        LOG4CXX_INFO(fileLogger, "domain-range x: " << dx);

        string dy = s1->settings().domain_size().y().get();
        domain_size[1] = std::stod(dy);
        LOG4CXX_INFO(fileLogger, "domain-range y: " << dy);

        string dz = s1->settings().domain_size().z().get();
        domain_size[2] = std::stod(dz);
        LOG4CXX_INFO(fileLogger, "domain-range z: " << dz);

        string rt = s1->settings().rcutoff().val().get();
        rcutoff = std::stod(rt);
        LOG4CXX_INFO(fileLogger, "cuttoff-radius: " << rt);

        string x1bc = s1->settings().x1_boundary_condition().val().get();
        if(x1bc == "outflow"){
            x1_boundary_condition = 1;
            LOG4CXX_INFO(fileLogger, "x1 boundary condition: outflow");
        }else if(x1bc=="reflecting") {
            x1_boundary_condition = 2;
            LOG4CXX_INFO(fileLogger, "x1 boundary condition: reflecting");
        }else if(x1bc=="periodic"){
            x1_boundary_condition = 3;
            LOG4CXX_INFO(fileLogger, "x1 boundary condition: periodic");
        }else{
            LOG4CXX_INFO(fileLogger, "" << x1bc << "" << " is not a valid boundary conditon");
        }

        string x2bc = s1->settings().x2_boundary_condition().val().get();
        if(x2bc == "outflow"){
            x2_boundary_condition = 1;
            LOG4CXX_INFO(fileLogger, "x2 boundary condition: outflow");
        }else if(x2bc=="reflecting") {
            x2_boundary_condition = 2;
            LOG4CXX_INFO(fileLogger, "x2 boundary condition: reflecting");
        }else if(x2bc=="periodic"){
            x2_boundary_condition = 3;
            LOG4CXX_INFO(fileLogger, "x2 boundary condition: periodic");
        }else{
            LOG4CXX_INFO(fileLogger, "" << x2bc << "" << " is not a valid boundary conditon");
        }

        string y1bc = s1->settings().y1_boundary_condition().val().get();
        if(y1bc == "outflow"){
            y1_boundary_condition = 1;
            LOG4CXX_INFO(fileLogger, "y1 boundary condition: outflow");
        }else if(y1bc=="reflecting") {
            y1_boundary_condition = 2;
            LOG4CXX_INFO(fileLogger, "y1 boundary condition: reflecting");
        }else if(y1bc=="periodic"){
            y1_boundary_condition = 3;
            LOG4CXX_INFO(fileLogger, "y1 boundary condition: periodic");
        }else{
            LOG4CXX_INFO(fileLogger, "" << y1bc << "" << " is not a valid boundary conditon");
        }

        string y2bc = s1->settings().y2_boundary_condition().val().get();
        if(y2bc == "outflow"){
            y2_boundary_condition = 1;
            LOG4CXX_INFO(fileLogger, "y2 boundary condition: outflow");
        }else if(y2bc=="reflecting") {
            y2_boundary_condition = 2;
            LOG4CXX_INFO(fileLogger, "y2 boundary condition: reflecting");
        }else if(y2bc=="periodic"){
            y2_boundary_condition = 3;
            LOG4CXX_INFO(fileLogger, "y2 boundary condition: periodic");
        }else{
            LOG4CXX_INFO(fileLogger, "" << y2bc << "" << " is not a valid boundary conditon");
        }

        string z1bc = s1->settings().z1_boundary_condition().val().get();
        if(z1bc == "outflow"){
            z1_boundary_condition = 1;
            LOG4CXX_INFO(fileLogger, "z1 boundary condition: outflow");
        }else if(z1bc=="reflecting") {
            z1_boundary_condition = 2;
            LOG4CXX_INFO(fileLogger, "z1 boundary condition: reflecting");
        }else if(z1bc=="periodic"){
            z1_boundary_condition = 3;
            LOG4CXX_INFO(fileLogger, "z1 boundary condition: periodic");
        }else{
            LOG4CXX_INFO(fileLogger, "" << z1bc << "" << " is not a valid boundary conditon");
        }

        string z2bc = s1->settings().z2_boundary_condition().val().get();
        if(z2bc == "outflow"){
            z2_boundary_condition = 1;
            LOG4CXX_INFO(fileLogger, "z2 boundary condition: outflow");
        }else if(z2bc=="reflecting") {
            z2_boundary_condition = 2;
            LOG4CXX_INFO(fileLogger, "z2 boundary condition: reflecting");
        }else if(z2bc=="periodic"){
            z2_boundary_condition = 3;
            LOG4CXX_INFO(fileLogger, "z2 boundary condition: periodic");
        }else{
            LOG4CXX_INFO(fileLogger, "" << z2bc << "" << " is not a valid boundary conditon");
        }

        string bm = s1->settings().brownian_motion().bol().get();
        if(bm == "true"){
            brownian_motion = true;
        }else if(bm == "false"){
            brownian_motion = false;
        }
        LOG4CXX_INFO(fileLogger, "brownian motion: " << bm);

        string it = s1->settings().initial_temperature().val().get();
        initial_temperature = std::stod(it);
        LOG4CXX_INFO(fileLogger, "initial temperature: " << it);

        string nt = s1->settings().n_thermostat().t().get();
        nthermostat = std::stod(nt);
        LOG4CXX_INFO(fileLogger, "nthermostat: " << nt);

        string tt = s1->settings().target_temperature().val().get();
        target_temperature = std::stod(tt);
        LOG4CXX_INFO(fileLogger, "target temperature: " << tt);

        string td = s1->settings().temperature_difference().val().get();
        temperature_difference = std::stod(td);
        LOG4CXX_INFO(fileLogger, "temperature difference: " << td);

        string g = s1->settings().gravitation().val().get();
        gravitation = std::stod(g);
        LOG4CXX_INFO(fileLogger, "gravitational force: " << g);

        string r00 = s1->settings().r0().val().get();
        r0 = std::stod(r00);
        LOG4CXX_INFO(fileLogger, "r0: " << r00);

        string k00 = s1->settings().k().val().get();
        k = std::stod(k00);
        LOG4CXX_INFO(fileLogger, "k: " << k00);

        string fz = s1->settings().fz_up().val().get();
        fz_up = std::stod(fz);
        LOG4CXX_INFO(fileLogger, "fz_up: " << fz);

        string pm = s1->settings().parallelisation_method().val().get();
        if(pm == "none"){
            parallelisation_method = 0;
            LOG4CXX_INFO(fileLogger, "parallelisation method: none");
        }else if(pm == "CoolMuc2"){
            parallelisation_method = 1;
            LOG4CXX_INFO(fileLogger, "parallelisation method: CoolMuc2");
        }else if(pm == "CoolMuc3"){
            parallelisation_method = 2;
            LOG4CXX_INFO(fileLogger, "parallelisation method: CoolMuc3");
        }else{
            LOG4CXX_INFO(fileLogger, "parallelisation method: none");
        }

        string fcm = s1->settings().force_calculation_method().val().get();
        if(fcm == "LJ"){
            force_calculation_method = 0;
            LOG4CXX_INFO(fileLogger, "force-calculation-method: LJ");
        }else if(fcm == "Smoothed LJ"){
            force_calculation_method = 1;
            LOG4CXX_INFO(fileLogger, "force-calculation-method: Smoothed LJ");
        }else if(fcm == "gravitation"){
            force_calculation_method = 1;
            LOG4CXX_INFO(fileLogger, "force-calculation-method: gravitation");
        }else{
            LOG4CXX_INFO(fileLogger, fcm << "is not a valid force-calculation-method, using standard");
        }

        string rlv = s1->settings().rl().val().get();
        rl = std::stod(rlv);
        LOG4CXX_INFO(fileLogger, "rl: " << rlv);

        //checkpoint
        string wc = s1->checkpoint()->write_checkpoint().bol().get();
        if(wc == "true"){
            writeCheckpoint = true;
            LOG4CXX_INFO(fileLogger, "writeCheckpoint: " << wc);
        }else if(wc == "false"){
            writeCheckpoint = false;
            LOG4CXX_INFO(fileLogger, "writeCheckpoint: " << wc);
        }else{
            LOG4CXX_INFO(fileLogger, "writeCheckpoint value has to be true/false. (" << wc << ") is not valid");
        }

        string rc = s1->checkpoint()->read_checkpoint().bol().get();
        if(rc == "true"){
            readCheckpoint = true;
            LOG4CXX_INFO(fileLogger, "readCheckpoint: " << rc);
        }else if(rc == "false"){
            readCheckpoint = false;
            LOG4CXX_INFO(fileLogger, "readCheckpoint: " << rc);
        }else {
            LOG4CXX_INFO(fileLogger, "readCheckpoint value has to be true/false. (" << rc << ") is not valid");
        }

        string wct = s1->checkpoint()->write_checkpoint_time().t().get();
        writeCheckpointTime = std::stod(wct);
        LOG4CXX_INFO(fileLogger, "write checkpoint time: " << wct);


        //spheres
        //double vBM = factor * sqrt(8 / 3.14);
        ///@brief take all spheres and initialize for each a particle generator with their characteristics
        for(shapesType::sphere_iterator s (s1->sphere().begin());
            s != s1->sphere ().end(); s++){
            x[0] = stod(s->centre().x().get())-domain_size[0]/2;
            x[1] = stod(s->centre().y().get())-domain_size[1]/2;
            x[2] = stod(s->centre().z().get())-domain_size[2]/2;
            r = stod(s->radius().val().get());
            h = stod(s->h().val().get());
            m = stod(s->mass().val().get());
            t = stod(s->type().val().get());
            e=stod(s->epsilon().val().get());
            si=stod(s->sigma().val().get());
            v[0] = stod(s->velocity().x().get());
            v[1] = stod(s->velocity().y().get());
            v[2] = stod(s->velocity().z().get());

            particleGenerator p(x, r, h, m, v, t, e, si);
            a.push_back(p);
        }
        ///@brief take all cuboids and initialize for each a particle generator with their characteristics
        for(shapesType::cuboid_const_iterator c (s1->cuboid().begin());
            c != s1->cuboid().end(); c++){

            x[0] = stod(c->right_top_point().x().get())-domain_size[0]/2;
            x[1] = stod(c->right_top_point().y().get())-domain_size[1]/2;
            x[2] = stod(c->right_top_point().z().get())-domain_size[2]/2;
            N[0] = stoi(c->side_lengths().x().get());
            N[1] = stoi(c->side_lengths().y().get());
            N[2] = stoi(c->side_lengths().z().get());
            h = stod(c->h().val().get());
            m = stod(c->mass().val().get());
            t = stod(c->type().val().get());
            e=stod(c->epsilon().val().get());
            si=stod(c->sigma().val().get());
            v[0] = stod(c->velocity().x().get());
            v[1] = stod(c->velocity().y().get());
            v[2] = stod(c->velocity().z().get());
            particleGenerator p(x, N, h, m, v, t, e, si);
            a.push_back(p);
        }
        return SettingsContainer (end_time, delta_t, factor, domain_size, rcutoff,
                                  brownian_motion, initial_temperature,
                                  nthermostat, target_temperature, temperature_difference, gravitation,
                                  writeCheckpoint, readCheckpoint, writeCheckpointTime,
                                  x1_boundary_condition, x2_boundary_condition, y1_boundary_condition,
                                  y2_boundary_condition, z1_boundary_condition, z2_boundary_condition,
                                  r0, k, fz_up, parallelisation_method, force_calculation_method, rl);

    }catch(const xml_schema::exception& e){
        LOG4CXX_ERROR(fileLogger, "Error reading file: "<< e);
        exit(-1);
    }
}



