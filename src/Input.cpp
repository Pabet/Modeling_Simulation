
#include <cstdlib>
#include <cstring>
#include <iostream>
#include "Input.h"
#include "FileReader.h"
#include <fstream>

#ifdef logOn
#include "log4cxx/logger.h"
#include "log4cxx/helpers/exception.h"
using namespace log4cxx;
using namespace log4cxx::helpers;
LoggerPtr inputLogger(Logger::getLogger("MolSim.Input"));
#else
#define LOG4CXX_TRACE(logger, expression)
#define LOG4CXX_DEBUG(logger, expression)
#define LOG4CXX_INFO(logger, expression)
#define LOG4CXX_WARN(logger, expression)
#define LOG4CXX_ERROR(logger, expression)
#define LOG4CXX_FATAL(logger, expression)
#endif

using namespace std;
Input::Input() {}

double end_timei ;
double delta_ti ;

std::list<Particle> Input::InterpretInput(int n, char* param[])
{
    ///@brief Input interpreter for txt files for the first and second task
    end_timei=1000;
    delta_ti=0.014;
    std::list<Particle> readparticles;
    ifstream f;
    f.open(param[2]);
    bool validfile=f.is_open();
    f.close();
    LOG4CXX_INFO(inputLogger, "Hello from MolSim for PSE!");
    string inputstr="eingabe-sonne.txt";
    char* input= &inputstr[0u];
    LOG4CXX_INFO(inputLogger, "The first argument of the program is the input file");
    LOG4CXX_INFO(inputLogger, "The second argument of the program is the end time");
    LOG4CXX_INFO(inputLogger, "The third argument of the program is the delta time");
    LOG4CXX_INFO(inputLogger, "In case of invalid input the default input will be taken");
    LOG4CXX_INFO(inputLogger, "The default input file is: "<< input);
    LOG4CXX_INFO(inputLogger, "The default end time is "<<end_timei<<" and the default delta time is "<<delta_ti);

    if (n != 5) {
        LOG4CXX_ERROR(inputLogger, "Errounous programme call! The default input will be taken");
    }
    else{
        if(validfile)
            input=param[2];
        else{LOG4CXX_ERROR(inputLogger, "The input file is invalid, default will be used");}
        LOG4CXX_INFO(inputLogger, "The input file is "<<input);
        double tempend_time=atof(param[3]);
        double tempdelta_t=atof(param[4]);
        if(tempend_time<=0.0)
        {LOG4CXX_ERROR(inputLogger, "Endtime is invalid, default will be used");}
        else
        {end_timei=tempend_time;}
        LOG4CXX_INFO(inputLogger, "End time is "<<end_timei);

        if(tempdelta_t<=0.0)
        {  LOG4CXX_ERROR(inputLogger, "Delta t is invalid, default will be used");}
        else
        {delta_ti=tempdelta_t;}
        LOG4CXX_INFO(inputLogger, "Delta t is "<<delta_ti);
    }

    FileReader fileReader;
    fileReader.readFile(readparticles, input);
    return readparticles;
}
