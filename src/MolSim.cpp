#include "outputWriter/VTKWriter.h"
#include "outputWriter/XYZWriter.h"
#include "Input.h"
#include "FileReader.h"
#include "outputWriter/pvdWriter.h"
//#include "MaxwellBoltzmannDistribution.h"
#include "particleGenerator.h"
#include "SettingsContainer.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <list>
#include <time.h>
#include <string>
#include <iomanip>

#include "Calculator.h"
#include "particleGenerator2.h"
#include "thermostat.h"
#include "outputWriter/checkpoint.h"
#include "Mesh.h"


#ifdef logOn
#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"
#include "log4cxx/helpers/exception.h"
#include <log4cxx/logger.h>
#include <log4cxx/xml/domconfigurator.h>

using namespace log4cxx;
using namespace log4cxx::xml;
using namespace log4cxx::helpers;
#else
#define LOG4CXX_TRACE(logger, expression)
#define LOG4CXX_DEBUG(logger, expression)
#define LOG4CXX_INFO(logger, expression)
#define LOG4CXX_WARN(logger, expression)
#define LOG4CXX_ERROR(logger, expression)
#define LOG4CXX_FATAL(logger, expression)
#endif

using namespace std;

///@image html graph1.png
///@image html graph2.png
/**** forward declaration of the calculation functions ****/

/**
 * calculate the force for all particles
 */
//void calculateF();

/**
 * calculate the position for all particles
 */
//void calculateX(std::list<Particle> &particles, double delta_t);

/**
 * calculate the position for all particles
 */
//void calculateV();

/**
 * plot the particles to a xyz-file
 */
//void calculateForce();
/**
 * calculate the new force for all particles
 */

 void plotParticles(int iteration);

double start_time = 0;
double end_time = 10;
double delta_t = 1;
Calculator calculator;
std::list<particleGenerator> a;
std::list<utils::Vector<double,3>> b;
double epsilon;
double sigma;
double factor;
std::list<Particle> particles;
utils::Vector<double, 3> domain_size;
double rcutoff;
int x1border;
int x2border;
int y1border;
int y2border;
int  z1border;
int z2border;
bool membrane=false;
double r0;
    double k;
    double fz_up;
    double rl;

double gravity;
double t_init=0;
double t_target=0;
double delta_temperature=0;
int n_thermo=0;
bool t_initialize=false;
double writeCheckpointTime=0;
bool wCheckpoint=false;
bool rCheckpoint=false;
int parallelisationMethod = 0;


int main(int argc, char *argsv[]) {
    ///@brief In order to make easy to change I/O functions , I/O will be put in a separate function

    // the forces are needed to calculate x, but are not given in the input file.
#ifdef logOn
    BasicConfigurator::configure();
    LoggerPtr logger(Logger::getLogger("MolSim"));
    LoggerPtr iterationLogger(Logger::getLogger("MolSim.Iteration"));
    DOMConfigurator::configure("Log4cxxConfig.xml");
#endif
    clock_t start = clock();
    string outputanim = "pvdout.pvd";
    outputWriter::pvdWriter collectionWriter;
    std::list<string> record;
    FileReader fileReader;
    SettingsContainer s = fileReader.readXMLFile(a, argsv[1]);
    domain_size = s.getDomain_size();
    rcutoff = s.getRcutoff();
    n_thermo = s.getNthermostat();
    end_time = s.ent_time;
    delta_t = s.delta_t;
    t_init = s.getInitial_temperature();
    t_target = s.getTarget_temperature();
    delta_temperature = s.getTemperature_difference();
    rCheckpoint = s.getReadCheckpoint();
    wCheckpoint = s.getWriteCheckpoint();
    writeCheckpointTime = s.getWriteCheckpointTime();
    parallelisationMethod = s.parallelisation_method;

    x1border = s.x1_boundary_condition;
    x2border = s.x2_boundary_condition;
    y1border = s.y1_boundary_condition;
    y2border = s.y2_boundary_condition;
    z1border = s.z1_boundary_condition;
    z2border = s.z2_boundary_condition;
r0=s.r0;
k=s.k;
fz_up=s.fz_up;
rl=s.rl;
    LOG4CXX_INFO(logger, "end time: " << end_time << "\tnext: " << delta_t << "\t" << endl);
    list<particleGenerator>::iterator iterator = a.begin();
    std::vector<Particle> particles2;
    if (n_thermo != 0) t_initialize = true;
    while (iterator != a.end()) {
        particleGenerator &f = *iterator;
        if (f.isCircle == true)
            f.createCircle(particles, 0.1, 2, t_init, t_initialize);
        else
            f.createGrid(particles, 0.1, t_init, 2, t_initialize,membrane,particles2);
        iterator++;
        LOG4CXX_INFO(logger, "size: " << particles.size() << endl);
    }


    list<utils::Vector<double, 3>>::iterator iteratora = b.begin();
    utils::Vector<double, 3> &a43 = *iteratora;
    epsilon = a43[0];
    sigma = a43[1];
    gravity = s.getGravitation();
    Mesh CellMesh
if (!membrane)
    CellMesh = Mesh(t_initialize,t_init,rcutoff, domain_size, particles, delta_t, gravity, x1border, x2border, y1border, y2border,
                         z1border, z2border, parallelisationMethod,rl);
else CellMesh = Mesh(t_initialize,t_init,rcutoff, domain_size, particles2, delta_t, gravity, x1border, x2border, y1border, y2border,
                         z1border, z2border,k,r0,fz_up,rl, parallelisationMethod);
    CellMesh.reassignMesh();
    //particleGenerator2 partgen(voutflow,vperiodic,vreflecting,houtflow,houtflow,hreflecting, domain_size, rcutoff, 1);


    //int retval = partgen.createMesh(particles, sigma, epsilon,true,gravity);


    double current_time = start_time;
    int iteration = 0;
    rCheckpoint = false;
    if (rCheckpoint) {

        readCheckpoint(particles, argsv[2], current_time);
    }
    //outputWriter::VTKWriter *a = new outputWriter::VTKWriter();
    //a->initializeOutput(particles.size());//Gehort zu MBD
bool rdf;
bool update;
    // for this loop, we assume: current x, current f and current v are known
    while (current_time < end_time) {
        // calculate new x
        if (rdf%1000==0) rdf=true;
        if (current_time<150) update=true;
        CellMesh.calculateX(rdf);

     //   double duration = (clock()-start)/(double) CLOCKS_PER_SEC;
     //   cout<<"runtime: "<<duration<<endl;

        CellMesh.halloCellsDecisions();
        // calculate new f
        CellMesh.calculateF(update);
rdf=false;
update=false;
        iteration++;

        // calculate new v
        CellMesh.calculateV();
        //if (wCheckpoint && (abs(current_time - writeCheckpointTime) < 0.00001)) {
        //    writeCheckpoint(particles, argsv[2], current_time);
       // }
        if (n_thermo != 0 && iteration % n_thermo == 0) {

            if (t_init != t_target && (((delta_temperature * (iteration / n_thermo)) + t_init) <= t_target))

                CellMesh.Target_Temperature(t_init + iteration * (delta_temperature / n_thermo));
            else CellMesh.Target_Temperature(t_target);
        }
        int y1=CellMesh.ysize;
        int z1=CellMesh.zsize;
        if (iteration % 25== 0) {
            outputWriter::VTKWriter writer1;
            writer1.initializeOutput(particles.size());
            //int sum=0;
            for (int i = 0; i < CellMesh.xsize; i++) {
                for (int j = 0; j < CellMesh.ysize; j++) {
                    for (int l = 0; l < CellMesh.zsize; l++) {
                        for (int k = 0; k < CellMesh.Cells[i*y1*z1+j*z1+l].getElementNumber(); k++) {
                            Particle* p = CellMesh.Cells[i*y1*z1+j*z1+l].getParticleAt(k);
                            writer1.plotParticle(*p);
                        }
                    }
                }
            }
            stringstream strstr;
            strstr << "VTKout" << "_" << setfill('0') << setw(4) << iteration << ".vtu";
            string str = strstr.str();
            stringstream current;
            current << "<DataSet timestep=\"" << (current_time) << "\" part=\"001\" file=\"" << "VTKout" << "_"
                    << setfill('0') << setw(4) << iteration << ".vtu" << "\"/>" << endl;
            record.push_back(current.str());
            writer1.writeFile("VTKout", iteration);
            //cout<<"written "<<iteration<<endl;
            //plotParticles(iteration);
        }
        LOG4CXX_DEBUG(iterationLogger, "" << iteration);
        current_time += delta_t;
    }
    collectionWriter.writeFile(outputanim, record);

    LOG4CXX_DEBUG(logger, "output written. Terminating...");

    clock_t end = clock();

    double elapsed = (double) (end - start) / CLOCKS_PER_SEC;
    ofstream fid;
    fid.open("CUTOFF8000;2secs.txt");
    fid << elapsed << endl;
    return 0;
}
void plotParticles(int iteration) {

  string out_name("MD_vtk");
  outputWriter::VTKWriter writer1;

  outputWriter::XYZWriter writer;
  writer.plotParticles(particles, out_name, iteration);
}
