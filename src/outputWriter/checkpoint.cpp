//
// Created by nastya on 18.12.18.
//

//
// Created by nastya on 17.12.18.
//
#include <list>
#include "checkpoint.h"
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <particleGenerator.h>

#ifdef logOn
#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"
#include "log4cxx/helpers/exception.h"
using namespace log4cxx;
using namespace log4cxx::helpers;
LoggerPtr checkPointLogger(Logger::getLogger("MolSim.Checkpoint"));
#endif
//FileReader::FileReader() {}

//FileReader::~FileReader() {}


using namespace std;

void writeCheckpoint (std::list<Particle> particles,const std::string &filename, double iteration){

    std::ofstream file;
    stringstream strstr;
    strstr << filename << "_file" << setfill('0') << setw(4) << iteration << ".txt";
    file.open(strstr.str().c_str());
    file << particles.size() << endl;

    
    file<<iteration<<endl;
    list<Particle>::iterator iterator = particles.begin();
    while (iterator != particles.end()) {
         Particle &p = *iterator;
         utils::Vector<double, 3> x = p.getX();
        
         file.setf(ios_base::showpoint);

         for (int i = 0; i < 3; i++) {
             file << x[i] << " ";
         }
x=p.getF();
        for (int i = 0; i < 3; i++) {
            file << x[i] << " ";
        }
        x=p.getOldF();
        for (int i = 0; i < 3; i++) {
            file << x[i] << " ";
        }
        x=p.getV();
        for (int i = 0; i < 3; i++) {
            file << x[i] << " ";
        }
        double a=p.getM();
        file<<a<<" ";
         int d=p.getType();
         file<<d<<" ";
         a=p.geteps();
         file<<a<<" ";
         a=p.getsig();
         file<<a<<" ";
         file << endl;
         iterator++;
     }

    iterator = particles.begin();
   
    file.close();


}

void readCheckpoint (std::list<Particle> &particles,const std::string filename, double &iteration){
    int counter=0;
    double x[] = {0, 0, 0};
    double v[] = {0, 0, 0};
    double F[]={0,0,0};
    double m = 1;
    int h=0;
    double F_old[]={0,0,0};
double eps=0;
double sig=0;
    std::ifstream input_file(filename);
    string tmp_string;
    if (input_file.is_open()) {

        getline(input_file, tmp_string);
  
        istringstream datastream(tmp_string);
        datastream>> counter;
     
        getline(input_file, tmp_string);
    
       
   
        istringstream datastream5(tmp_string);
        datastream5>>iteration;
        getline(input_file, tmp_string);

        for (int i = 0; i < counter; i++) {
            istringstream datastream(tmp_string);
     
            for (int j = 0; j < 3; j++) {
                datastream >> x[j];
            }
            for (int j = 0; j < 3; j++) {
                datastream >> F[j];
            }
            for (int j=0; j<3;j++){
                datastream>>F_old[j];}

            for (int j=0; j<3;j++){
                datastream>>v[j];}

         
datastream>>m;
            datastream>>h;
            datastream>>eps;
            datastream>>sig;
               if (datastream.eof()) {
          cout
                  << "Error reading file: eof reached unexpectedly reading from line "
                 << i << endl;
          exit(-1);
     }
            Particle p(x,v,m,h,eps,sig);
            p.getF().operator=(F);
            p.getOldF().operator=(F_old);
            particles.push_back(p);



            getline(input_file, tmp_string);
       
        }


   




    } else {
        std::cout << "Error: could not open file " << filename << std::endl;
        exit(-1);
    }
}