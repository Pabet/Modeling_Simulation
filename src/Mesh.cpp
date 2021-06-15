#include "Mesh.h"
#include "particleContainer.h"
#include "utils/Vector.h"
#include "Particle.h"
#include "Calculator.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef logOn
#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"
#include "log4cxx/helpers/exception.h"
using namespace log4cxx;
using namespace log4cxx::helpers;
LoggerPtr threadLogger(Logger::getLogger("MolSim.Threads"));
#else
#define LOG4CXX_TRACE(logger, expression)
#define LOG4CXX_DEBUG(logger, expression)
#define LOG4CXX_INFO(logger, expression)
#define LOG4CXX_WARN(logger, expression)
#define LOG4CXX_ERROR(logger, expression)
#define LOG4CXX_FATAL(logger, expression)
#endif

Mesh::Mesh(bool initTemp,double t_Target,double rcutoff,utils::Vector<double, 3> domain_size,std::list<Particle>&tempparticles,double d_t,double gravitation,int x1bound,int x2bound,int y1bound,int y2bound,int z1bound,int z2bound, int pM,double r2)
{membrane=false;
 rl=r2;
    inittemp=initTemp;
    t_target=t_Target;
    x1boundary=x1bound;
    x2boundary=x2bound;
    y1boundary=y1bound;
    y2boundary=y2bound;
    z1boundary=z1bound;
    z2boundary=z2bound;
    s=Calculator();
    delta_t=d_t;
    gravitationa=gravitation;
    zdim=true;
    cutoff=rcutoff;
    xsize=ceil(domain_size.operator[](0)/cutoff)+2;
    ysize=ceil(domain_size.operator[](1)/cutoff)+2;
    zsize=ceil(domain_size.operator[](2)/cutoff)+2;
    if(domain_size.operator[](2)==0) {
        zsize = 3;
        zdim = false;
        z1boundary=0;
        z2boundary=0;
    }
    this->domain_size=utils::Vector<double, 3>(domain_size);
    Cells=new particleContainer[xsize*ysize*zsize];
    particles=tempparticles;
    parallelisationMethod=pM;
}
Mesh::Mesh(bool initTemp,double t_Target,double rcutoff,utils::Vector<double, 3> domain_size,std::vector<Particle>&tempparticles,double d_t,double gravitation,int x1bound,int x2bound,int y1bound,int y2bound,int z1bound,int z2bound,double k1, double r1, double F_z,double r2,int pM)
{k=k1;
r0=r1;
F_z_UP=F_z;
membrane=true;
rl=r2;
inittemp=initTemp;
    t_target=t_Target;
    x1boundary=x1bound;
    x2boundary=x2bound;
    y1boundary=y1bound;
    y2boundary=y2bound;
    z1boundary=z1bound;
    z2boundary=z2bound;
    s=Calculator();
    delta_t=d_t;
    gravitationa=gravitation;
    zdim=true;
        parallelisationMethod=pM;
    cutoff=rcutoff;
    xsize=ceil(domain_size.operator[](0)/cutoff)+2;
    ysize=ceil(domain_size.operator[](1)/cutoff)+2;
    zsize=ceil(domain_size.operator[](2)/cutoff)+2;
    if(domain_size.operator[](2)==0) {
        zsize = 3;
        zdim = false;
        z1boundary=0;
        z2boundary=0;
    }
    this->domain_size=utils::Vector<double, 3>(domain_size);
    Cells=new particleContainer[xsize*ysize*zsize];
    particles2=tempparticles;
     parallelisationMethod=pM;
}
double diffusion=0.0;
void Mesh::reassignMesh()
{
    std::list<Particle>::iterator iterator;
    iterator=particles.begin();
    while(iterator!=particles.end())
    {if (!membrane) {
        Particle &p1=*iterator;
        Particle* p2=new Particle(p1);
        double tx=double(xsize)/2;
        double ty=double(ysize)/2;
        double tz=double(zsize)/2;
        int xposition= floor(p1.getX().operator[](0)/cutoff+tx);
        int yposition= floor(p1.getX().operator[](1)/cutoff+ty);
        int zposition= floor(p1.getX().operator[](2)/cutoff+tz);
        if(!zdim)
            zposition=1;
        if((x1boundary==1&&xposition<=0)||(x2boundary==1&&xposition>=xsize-1)||(y1boundary==1&&yposition<=0)||(y2boundary==1&&yposition>=ysize-1)||(zdim&&((z1boundary==1&&zposition<=0)||(z2boundary==1&&zposition>=zsize-1))))
        {
            iterator=particles.erase(iterator);
        }
        else {
            if (xposition >= 0 && yposition >= 0 && xposition < xsize && yposition < ysize &&
                (!zdim || (zposition >= 0 && zposition < zsize))) {
                Cells[xposition*ysize*zsize+yposition*zsize+zposition].addParticle(p2);
            }
            iterator++;
        }
    }}
    else {  std::vector<Particle>::iterator iterator;
        iterator=particles2.begin();
        while(iterator!=particles2.end())
        {
            Particle &p1=*iterator;
            Particle* p2=new Particle(p1);
            double tx=double(xsize)/2;
            double ty=double(ysize)/2;
            double tz=double(zsize)/2;
            int xposition= floor(p1.getX().operator[](0)/cutoff+tx);
            int yposition= floor(p1.getX().operator[](1)/cutoff+ty);
            int zposition= floor(p1.getX().operator[](2)/cutoff+tz);
            if(!zdim)
                zposition=1;
            if((x1boundary==1&&xposition<=0)||(x2boundary==1&&xposition>=xsize-1)||(y1boundary==1&&yposition<=0)||(y2boundary==1&&yposition>=ysize-1)||(zdim&&((z1boundary==1&&zposition<=0)||(z2boundary==1&&zposition>=zsize-1))))
            {
                iterator=particles2.erase(iterator);
            }
            else {
                if (xposition >= 0 && yposition >= 0 && xposition < xsize && yposition < ysize &&
                    (!zdim || (zposition >= 0 && zposition < zsize))) {
                    Cells[xposition*ysize*zsize+yposition*zsize+zposition].addParticle(p2);
                }
                iterator++;
            }
        }}
    if(inittemp)
    {
        this->Target_Temperature(t_target);
    }
}
void Mesh::calculateF(bool timeF)
{
    particleContainer* tempx= new particleContainer[ysize*zsize];
    particleContainer* tempy= new particleContainer[xsize*zsize];
    particleContainer* tempz= new particleContainer[xsize*ysize];

    for(int i1=1;i1<xsize-1;i1++)
    {
        for(int j1=1;j1<ysize-1;j1++)
        {
            for(int l1=1;l1<zsize-1;l1++)
            {
                for(int k1=0;k1<Cells[i1*ysize*zsize+j1*zsize+l1].getElementNumber();k1++)
                {
                    utils::Vector<double ,3> summedForce=utils::Vector<double ,3>(0.0);
                    Particle &p1=*Cells[i1*ysize*zsize+j1*zsize+l1].getParticleAt(k1);
                    double distance=(pow(2,1/6)*p1.getsig())+1;
                    double d_new []={0,0,0} ;
                    utils::Vector<double ,3> distanceVector=utils::Vector<double ,3>(d_new);///@brief if distance is not changed no changes will be made to force
                    if(x1boundary==2&&i1==1) ///@brief calculation of distance to a pseudo molecule for a reflective border
                    {
                        distance=2*(p1.getX().operator[](0)+(domain_size.operator[](0)/2));
                        d_new[0] =-distance;
                        distanceVector=utils::Vector<double ,3>(d_new);
                    }
                    if(x2boundary==2&&i1==(xsize-2))
                    {
                        distance=2*(p1.getX().operator[](0)-(domain_size.operator[](0)/2));
                        d_new [0]=-distance;
                        distanceVector=utils::Vector<double ,3>(d_new);
                    }
                    if(x1boundary==2||x2boundary==2) {
                        if (distance < (pow(2, 1.0 / 6.0) * p1.getsig())) {
                            double sigmaMean = p1.getsig() * p1.getsig();
                            double epsilonMean = p1.geteps() + p1.geteps();
                            double potentialplus = ((12 * epsilonMean * pow(sigmaMean, 3)) / pow((distance), 8));
                            double potentialminus = -((pow(sigmaMean, 6) * 24 * epsilonMean) / pow((distance), 14));
                            utils::Vector<double, 3> force = (distanceVector.operator*(
                                    potentialplus)).operator+(distanceVector.operator*(potentialminus));
                            summedForce.operator=(summedForce.operator+(force));
                        }
                    }
                    d_new [0]=0 ;
                    if(y1boundary==2&&j1==1)
                    {
                        distance=2*(p1.getX().operator[](1)+(domain_size.operator[](1)/2));
                        d_new [1]=-distance;
                        distanceVector=utils::Vector<double ,3>(d_new);
                    }
                    if(y2boundary==2&&j1==(ysize-2))
                    {
                        distance=2*(p1.getX().operator[](1)-(domain_size.operator[](1)/2));
                        d_new [1]=-distance;
                        distanceVector=utils::Vector<double ,3>(d_new);
                    }
                    if(y1boundary==2||y2boundary==2) {
                        if (distance < (pow(2, 1.0 / 6.0) * p1.getsig())) {
                            double sigmaMean = p1.getsig() * p1.getsig();
                            double epsilonMean = p1.geteps() + p1.geteps();
                            double potentialplus = ((12 * epsilonMean * pow(sigmaMean, 3)) / pow((distance), 8));
                            double potentialminus = -((pow(sigmaMean, 6) * 24 * epsilonMean) / pow((distance), 14));
                            utils::Vector<double, 3> force = (distanceVector.operator*(
                                    potentialplus)).operator+(distanceVector.operator*(potentialminus));
                            summedForce.operator=(summedForce.operator+(force));
                        }
                    }
                    d_new [1]=0 ;
                    if(z1boundary==2&&l1==1)
                    {
                        distance=2*(p1.getX().operator[](2)+(domain_size.operator[](2)/2));
                        d_new [2]=-distance;
                        distanceVector=utils::Vector<double ,3>(d_new);
                    }
                    if(z2boundary==2&&l1==(zsize-2))
                    {
                        distance=2*(p1.getX().operator[](2)-(domain_size.operator[](2)/2));
                        d_new [2]=-distance;
                        distanceVector=utils::Vector<double ,3>(d_new);
                    }
                    if(z1boundary==2||z2boundary==2) {
                        if (distance < (pow(2, 1.0 / 6.0) * p1.getsig())) {
                            double sigmaMean = p1.getsig() * p1.getsig();
                            double epsilonMean = p1.geteps() + p1.geteps();
                            double potentialplus = ((12 * epsilonMean * pow(sigmaMean, 3)) / pow((distance), 8));
                            double potentialminus = -((pow(sigmaMean, 6) * 24 * epsilonMean) / pow((distance), 14));
                            utils::Vector<double, 3> force = (distanceVector.operator*(
                                    potentialplus)).operator+(distanceVector.operator*(potentialminus));
                            summedForce.operator=(summedForce.operator+(force));
                        }
                    }
                    for (int i2 = i1 - 1; i2 <= i1 ; i2++)
                    {
                        bool periodicborderx1=false;
                        bool periodicborderx2=false;
                        if(i2==0&&x1boundary==3)///@brief taking the neighboring periodic cell in consideration
                        {
                            i2=xsize-2;
                            periodicborderx1=true;
                        }
                        if(i2==xsize-1&&x2boundary==3)
                        {
                            i2=1;
                            periodicborderx2=true;
                        }
                        if(i2>0&&i2<xsize-1) {
                            for (int j2 = j1 - 1; j2 <= j1 + 1 ; j2++) {
                                if (i2 != i1 || j2 != j1 + 1) {
                                    bool periodicbordery1 = false;
                                    bool periodicbordery2 = false;
                                    if (j2 == 0 && y1boundary == 3) {
                                        j2=ysize-2;
                                        periodicbordery1 = true;
                                    }
                                    if (j2 == ysize - 1 && y2boundary == 3) {
                                        j2 = 1;
                                        periodicbordery2 = true;
                                    }
                                    if (j2 > 0 && j2 < xsize - 1) {
                                        for (int l2 = l1 - 1; l2 <= l1; l2++) {
                                            if(i1!=i2||j1!=i2||l2!=l1+1) {
                                                bool periodicborderz1 = false;
                                                bool periodicborderz2 = false;
                                                if (l2 == 0 && z1boundary == 3) {
                                                    l2 = zsize - 2;
                                                    periodicborderz1 = true;
                                                }
                                                if (l2 == zsize - 1 && z2boundary == 3) {
                                                    l2 = 1;
                                                    periodicborderz2 = true;
                                                }
                                                if (l2 > 0 && l2 < zsize - 1) {
                                                    for (int k2 = 0; k2 < Cells[i2*ysize*zsize+j2*zsize+l2].getElementNumber(); k2++) {
                                                        Particle &p2 = *Cells[i2*ysize*zsize+j2*zsize+l2].getParticleAt(k2);
                                                        if (!p1.operator==(p2)) {
                                                            utils::Vector<double, 3> position1 = utils::Vector<double, 3>(
                                                                    p1.getX());
                                                            utils::Vector<double, 3> position2 = utils::Vector<double, 3>(
                                                                    p2.getX());
                                                            utils::Vector<double, 3> distanceVector = utils::Vector<double, 3>(
                                                                    position2.operator-(position1));
                                                            if (periodicborderx1 ||
                                                                periodicborderx2)///@brief adjustment to distance of periodic boundaries
                                                            {
                                                                double d_new[] = {domain_size.operator[](0), 0, 0};
                                                                if (distanceVector.operator[](0) > 0) {
                                                                    distanceVector.operator=(utils::Vector<double, 3>(
                                                                            distanceVector.operator-(
                                                                                    utils::Vector<double, 3>(d_new))));
                                                                } else {
                                                                    distanceVector.operator=(utils::Vector<double, 3>(
                                                                            distanceVector.operator+(
                                                                                    utils::Vector<double, 3>(d_new))));
                                                                }
                                                            }
                                                            if (periodicbordery1 || periodicbordery2) {
                                                                double d_new[] = {0, domain_size.operator[](1), 0};
                                                                if (distanceVector.operator[](1) > 0) {
                                                                    distanceVector.operator=(utils::Vector<double, 3>(
                                                                            distanceVector.operator-(
                                                                                    utils::Vector<double, 3>(d_new))));
                                                                } else {
                                                                    distanceVector.operator=(utils::Vector<double, 3>(
                                                                            distanceVector.operator+(
                                                                                    utils::Vector<double, 3>(d_new))));
                                                                }
                                                            }
                                                            if (periodicborderz1 || periodicborderz2) {
                                                                double d_new[] = {0, 0, domain_size.operator[](2)};
                                                                if (distanceVector.operator[](2) > 0) {
                                                                    distanceVector.operator=(utils::Vector<double, 3>(
                                                                            distanceVector.operator-(
                                                                                    utils::Vector<double, 3>(d_new))));
                                                                } else {
                                                                    distanceVector.operator=(utils::Vector<double, 3>(
                                                                            distanceVector.operator+(
                                                                                    utils::Vector<double, 3>(d_new))));
                                                                }
                                                            }

                                                            double distance = distanceVector.L2Norm();

                                                            if (distance <= cutoff) {
                                                                utils::Vector<double,3>force(0.0);
                                                                if (membrane){

                                                                    Particle p3;
                                                                    if (p1.getneighbors()[0]!=NULL)
                                                                    {p3=*p1.getneighbors()[0];
                                                                    force.operator=((p3.getX().operator-(p1.getX())).operator*(k*((p1.getX().operator-(p3.getX())).L2Norm()-r0)));}
                                                                    if (p1.getneighbors()[1]!=NULL){
                                                                        p3=*p1.getneighbors()[1];
                                                                        force.operator=(force.operator+((p3.getX().operator-(p1.getX())).operator*(k*((p1.getX().operator-(p3.getX())).L2Norm()-r0))));
                                                                    }
if (p1.getneighbors()[2]!=NULL) {
    p3=*p1.getneighbors()[2];
    force.operator=(force.operator+((p3.getX().operator-(p1.getX())).operator*(k*((p1.getX().operator-(p3.getX())).L2Norm()-r0))));
}
if (p1.getneighbors()[3]!=NULL){p3=*p1.getneighbors()[3];
    force.operator=(force.operator+((p3.getX().operator-(p1.getX())).operator*(k*((p1.getX().operator-(p3.getX())).L2Norm()-r0))));}
if (p1.getDiagNeighbors()[0]!=NULL){
   p3=*p1.getDiagNeighbors()[0];
        force.operator=(force.operator+((p3.getX().operator-(p1.getX())).operator*(k*((p1.getX().operator-(p3.getX())).L2Norm()-(sqrt(2)*r0)))));
}
if (p1.getDiagNeighbors()[1]!=NULL){
    p3=*p1.getDiagNeighbors()[1];
    force.operator=(force.operator+((p3.getX().operator-(p1.getX())).operator*(k*((p1.getX().operator-(p3.getX())).L2Norm()-(sqrt(2)*r0)))));
}
if (p1.getDiagNeighbors()[2]!=NULL){
    p3=*p1.getDiagNeighbors()[2];
    force.operator=(force.operator+((p3.getX().operator-(p1.getX())).operator*(k*((p1.getX().operator-(p3.getX())).L2Norm()-(sqrt(2)*r0)))));
}
if (p1.getDiagNeighbors()[3]!=NULL){
    p3=*p1.getDiagNeighbors()[3];
    force.operator=(force.operator+((p3.getX().operator-(p1.getX())).operator*(k*((p1.getX().operator-(p3.getX())).L2Norm()-(sqrt(2)*r0)))));
}

if (distance<(pow(2,1/6)*p1.getsig())) {
    double sigmaMean = p1.getsig() * p2.getsig();
    double epsilonMean = p1.geteps() + p2.geteps();
    double potentialplus = (
            (12 * epsilonMean * pow(sigmaMean, 3)) /
            pow(distance, 8));
    double potentialminus = -(
            (pow(sigmaMean, 6) * 24 * epsilonMean) /
            pow(distance, 14));

   force.operator=(force.operator+ ((distanceVector.operator*(
            potentialplus)).operator+(
            distanceVector.operator*(potentialminus))));

}                         } else {if (rl==0 || (distance<=rl)) {
                                                                double sigmaMean = p1.getsig() * p2.getsig();
                                                                double epsilonMean = p1.geteps() + p2.geteps();
                                                                double potentialplus = (
                                                                        (12 * epsilonMean * pow(sigmaMean, 3)) /
                                                                        pow(distance, 8));
                                                                double potentialminus = -(
                                                                        (pow(sigmaMean, 6) * 24 * epsilonMean) /
                                                                        pow(distance, 14));

                                                                 force.operator=( (distanceVector.operator*(
                                                                        potentialplus)).operator+(
                                                                        distanceVector.operator*(potentialminus)));}
                                                                        else {
                                                                            double a=pow (p1.getsig(),6);
                                                                            double b=pow(distance,6);
                                                                            double c=cutoff-distance;
                                                                            force.operator=((p2.getX().operator-(p1.getX())).operator*((24*a*p1.geteps()/(distance*b*pow(cutoff-rl,3)))*((c*c
                                                                                    *(b-2*a)*(cutoff-3*rl+2*distance))-distance*(distance-cutoff)*(distance-rl)*(b-a))));

                                                                        }
                                                                if(!(i2==i1&&j2==j1&&l2==l1)) {
                                                                    p2.getF().operator=(p2.getF().operator-(force));
                                                                    if(periodicborderx1)
                                                                    {   Particle* newx=new Particle[1];
                                                                        newx[0]=Particle(p2.getX(),force,0,0,0,0);
                                                                        tempx[j1*zsize+l1].addParticle(newx);};
                                                                    if(periodicbordery1)
                                                                       {   Particle* newy=new Particle[1];
                                                                        newy[0]=Particle(p2.getX(),force,0,0,0,0);
                                                                        tempy[i1*zsize+l1].addParticle(newy);}
                                                                    if(periodicborderz1)
                                                                    {   Particle* newz=new Particle[1];
                                                                        newz[0]=Particle(p2.getX(),force,0,0,0,0);
                                                                        tempz[i1*ysize+j1].addParticle(newz);}
                                                                }}
                                                                summedForce.operator=(summedForce.operator+(force));
                                                            }
                                                        }
                                                    }
                                                }
                                                if (periodicborderz1 == true) {
                                                    l2 = 0;
                                                }
                                                if (periodicborderz2 == true) {
                                                    l2 = zsize - 1;
                                                }
                                            }
                                        }
                                    }
                                    if (periodicbordery1 == true) {
                                        j2 = 0;
                                    }
                                    if (periodicbordery2 == true) {
                                        j2 = ysize - 1;
                                    }
                                }
                            }
                        }
                        if(periodicborderx1== true)
                        {
                            i2=0;
                        }
                        if(periodicborderx2==true)
                        {
                            i2=xsize-1;
                        }
                    }
                    p1.getOldF().operator=(p1.getF());
                    double f_new[]={0,(p1.getM()*gravitationa),0};
                    utils::Vector<double,3>ga=utils::Vector<double,3>(f_new);
                    summedForce.operator=(summedForce.operator+(ga));///@brief addition of gravitational acceleration
                    p1.getF().operator=(summedForce);
                    summedForce=utils::Vector<double ,3>(0.0);
                    if(x1boundary==3&&i1==xsize-2)
                    {
                        for(int k3=0;k3<tempx[j1*zsize+l1].getElementNumber();k3++) {
                            Particle f=*tempx[j1 * zsize + l1].getParticleAt(k3);
                            if(f.getX().operator==(p1.getX()))
                            p1.getF().operator=(p1.getF().operator-(f.getV()));
                        }
                    }
                    if(y1boundary==3&&j1==ysize-2)
                        {
                        for(int k3=0;k3<tempy[i1*zsize+l1].getElementNumber();k3++) {
                            Particle f=*tempy[i1*zsize+l1].getParticleAt(k3);
                            if(f.getX().operator==(p1.getX()))
                            p1.getF().operator=(p1.getF().operator-(f.getV()));
                        }
                    }
                    if(z1boundary==3&&l1==zsize-2)
                    {
                        for(int k3=0;k3<tempz[i1*zsize+j1].getElementNumber();k3++) {
                            Particle f=*tempz[i1*zsize+j1].getParticleAt(k3);
                            if(f.getX().operator==(p1.getX()))
                                p1.getF().operator=(p1.getF().operator-(f.getV()));
                        }
                    }
                }
            }
        }

    }
    delete []tempx;
    delete []tempy;
    delete []tempz;
    if (membrane && timeF) {
        utils::Vector<double,3> force(0.0);
        force[2]=F_z_UP;
        particles2[17*50+24].getF().operator=(particles2[17*50+24].getF().operator+(force));
        particles2[17*50+25].getF().operator=(particles2[17*50+25].getF().operator+(force));
        particles2[18*50+24].getF().operator=(particles2[18*50+24].getF().operator+(force));
        particles2[18*50+25].getF().operator=(particles2[18*50+25].getF().operator+(force));
    }
}

void Mesh::calculateX(bool rdf) {
    std::list<Particle*> donelist;
    if (rdf) diffusion=0.0;
#ifdef _OPENMP
#pragma omp parallel for collapse(3) shared(donelist,Cells)
#endif
    for (int i = 1; i < xsize-1; i++) {
        for (int j = 1; j < ysize-1; j++) {
            for (int l = 1; l < zsize-1; l++) {
                for (int k = 0; k < Cells[i*ysize*zsize+j*zsize+l].getElementNumber(); k++) {
                    bool done=false;
                    Particle* p1 = Cells[i*ysize*zsize+j*zsize+l].getParticleAt(k);
                    std::list<Particle*>::iterator iterator;
                    iterator=donelist.begin();
                    while(iterator!=donelist.end())
                    {
                        Particle *p2=*iterator;
                        if(p2->operator==(*p1))
                        {
                            done=true;
                        }
                        iterator++;
                    }
                    if(!done) {
                        s.calculateX1P(*p1, delta_t);
                        if (rdf){
                      double f[]= {     p1->getX().operator[](0)+p1->getPeriodicMovement()[0]*domain_size[0],
                                                        p1->getX().operator[](1)+p1->getPeriodicMovement()[1]*domain_size[1],
                                                                p1->getX().operator[](2)+p1->getPeriodicMovement()[2]*domain_size[2]};
                      utils::Vector<double,3>realX(f);
                      double b=(realX.operator-(p1->getOldX())).L2Norm();
                      diffusion=diffusion+b*b;
                      p1->getOldX().operator=(p1->getX());
                        }
                        double tx=double(xsize)/2;
                        double ty=double(ysize)/2;
                        double tz=double(zsize)/2;
                        int xposition= floor(p1->getX().operator[](0)/cutoff+tx);
                        int yposition= floor(p1->getX().operator[](1)/cutoff+ty);
                        int zposition= floor(p1->getX().operator[](2)/cutoff+tz);
                        if(!zdim)
                            zposition=1;
                        if(xposition!=i||yposition!=j||zposition!=l)
                        {
                            if(xposition>=0&&xposition<xsize&&yposition>=0&&yposition<ysize&&zposition>=0&&zposition<zsize)
                            Cells[xposition*ysize*zsize+yposition*zsize+zposition].addParticle(p1);
                            Cells[i*ysize*zsize+j*zsize+l].removeparticleAt(k);
                            k--;
                        }
                        donelist.push_back(p1);
                    }

                }
            }
        }
    }
    donelist.clear();
}
void Mesh::calculateV() {
    if(parallelisationMethod==1){
        //CoolMuc2 14x2 Threads
        omp_set_num_threads(28);
    }else if(parallelisationMethod==2){
        //CoolMuc3
        omp_set_num_threads(256);
    }else{
        //no parallelisation
        omp_set_num_threads(1);
    }
#ifdef _OPENMP
#pragma omp parallel for collapse(3)
#endif
    for (int i = 1; i < xsize-1; i++) {
        for (int j = 1; j < ysize-1; j++) {
            for (int l = 1; l < zsize-1; l++) {
                for (int k = 0; k < Cells[i*ysize*zsize+j*zsize+l].getElementNumber(); k++) {
                    Particle *p1=Cells[i*ysize*zsize+j*zsize+l].getParticleAt(k);
                    s.calculateV1P(*p1,delta_t);
                }
            }
        }
    }
}
void Mesh::Target_Temperature (double t_target){
    int particlenumber=0;
    double e_kin=0;
     for (int j = 0; j < ysize; j++) {
         for (int l = 0; l < zsize; l++) {
             for (int k = 0; k < Cells[j * zsize + l].getElementNumber(); k++) {
                 Particle p = *Cells[j * zsize + l].getParticleAt(k);
                 particlenumber++;
                 e_kin = e_kin +
                         0.5 * (p.getV()[0] * p.getV()[0] + p.getV()[1] * p.getV()[1] + p.getV()[2] * p.getV()[2]);
             }
         }
     }
     int dimensions=3;
     if(!zdim)
         dimensions=2;
    double temp=2*e_kin/(dimensions*particlenumber);
    double betha=sqrt(t_target/temp);
    for (int j = 0; j < ysize; j++) {
        for (int l = 0; l < zsize; l++) {
            for (int k = 0; k < Cells[j * zsize + l].getElementNumber(); k++) {
                Particle p = *Cells[j * zsize + l].getParticleAt(k);
                p.getV().operator=(p.getV().operator*(betha));
            }
        }
    }

}
double Mesh::getRDF(){return diffusion/(domain_size[0]*domain_size[1]*domain_size[2]);}
void Mesh::halloCellsDecisions() {
    if(parallelisationMethod==1){
        //CoolMuc2 14x2 Threads
        omp_set_num_threads(28);
    }else if(parallelisationMethod==2){
        //CoolMuc3
        omp_set_num_threads(256);
    }else{
        //no parallelisation
        omp_set_num_threads(1);
    }
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < ysize; j++) {
        for (int l = 0; l < zsize; l++) {
            for (int k = 0; k < Cells[j*zsize+l].getElementNumber(); k++) {
                Particle* p = Cells[j*zsize+l].getParticleAt(k);
                switch (x1boundary) {
                    case 1: { ///@brief remove particle if outflow
                        Cells[j*zsize+l].removeparticleAt(k);
                        k--;
                        break;
                    }
                    case 2: { ///@brief if boundary reflective particle should not go out
                        break;
                    }
                    case 3: {
                        ///@brief set particle on the other side if periodic
                        double x_new[] = {(p->getX().operator[](0) + domain_size.operator[](0)), p->getX().operator[](1),
                                          p->getX().operator[](2)};
                                          utils::Vector<int,3> a(0.0);
                        a[0]=1;
                        p->getPeriodicMovement().operator+(a);
                        p->getX().operator=(utils::Vector<double, 3>(x_new));
                        Cells[(xsize - 2)*ysize*zsize+j*zsize+l].addParticle(p);
                        Cells[j*zsize+l].removeparticleAt(k);
                        k--;
                        break;
                    }
                    default: {
                        break;
                    }
                }
            }
            for (int k = 0; k < Cells[(xsize - 1)*ysize*zsize+j*zsize+l].getElementNumber(); k++) {
                Particle* p = Cells[(xsize - 1)*ysize*zsize+j*zsize+l].getParticleAt(k);
                switch (x2boundary) {
                    case 1: {///@brief remove particle if outflow
                        Cells[(xsize - 1)*ysize*zsize+j*zsize+l].removeparticleAt(k);
                        k--;
                        break;
                    }
                    case 2: ///@brief if boundary reflective particle should not go out
                    {
                        break;
                    }
                    case 3: {
                        ///@brief set particle on the other side if periodic
                        double x_new[] = {(p->getX().operator[](0) - domain_size.operator[](0)), p->getX().operator[](1),
                                          p->getX().operator[](2)};
                                          utils::Vector<int,3> a(0.0);
                        a[0]=-1;
                        p->getPeriodicMovement().operator+(a);
                        p->getX().operator=(utils::Vector<double, 3>(x_new));
                        Cells[ysize*zsize+j*zsize+l].addParticle(p);
                        Cells[(xsize - 1)*ysize*zsize+j*zsize+l].removeparticleAt(k);
                        k--;
                        break;
                    }
                    default:
                        break;
                }
            }
        }
    }
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int i = 0; i < xsize; i++) {
        for (int l = 0; l < zsize; l++) {
            for (int k = 0; k < Cells[i*ysize*zsize+l].getElementNumber(); k++) {
                Particle* p = Cells[i*ysize*zsize+l].getParticleAt(k);
                switch (y1boundary) {
                    case 1: { ///@brief remove particle if outflow
                        Cells[i*ysize*zsize+l].removeparticleAt(k);
                        k--;
                        break;
                    }
                    case 2: ///@brief if boundary reflective particle should not go out
                    {
                        break;
                    }
                    case 3: {///@brief set particle on the other side if periodic
                        double x_new[] = {p->getX().operator[](0), (p->getX().operator[](1) + domain_size.operator[](1)),
                                          p->getX().operator[](2)};
                                           utils::Vector<int,3> a(0.0);
                        a[1]=1;
                        p->getPeriodicMovement().operator+(a);
                        p->getX().operator=(utils::Vector<double, 3>(x_new));
                        Cells[i*ysize*zsize+(ysize - 2)*zsize+l].addParticle(p);
                        Cells[i*ysize*zsize+l].removeparticleAt(k);
                        k--;
                        break;
                    }
                    default: {
                        break;
                    }
                }
            }
            for (int k = 0; k < Cells[i*ysize*zsize+(ysize - 1)*zsize+l].getElementNumber(); k++) {
                Particle* p = Cells[i*ysize*zsize+(ysize - 1)*zsize+l].getParticleAt(k);
                switch (y2boundary) {
                    case 1: ///@brief remove particle if outflow
                    {
                        Cells[i*ysize*zsize+(ysize - 1)*zsize+l].removeparticleAt(k);
                        k--;
                        break;
                    }
                    case 2: ///@brief if boundary reflective particle should not go out
                    {
                        break;
                    }
                    case 3: {///@brief set particle on the other side if periodic
                        double x_new[] = {p->getX().operator[](0), (p->getX().operator[](1) - domain_size.operator[](1)),
                                          p->getX().operator[](2)};
                                            utils::Vector<int,3> a(0.0);
                        a[1]=-1;
                        p->getPeriodicMovement().operator+(a);
                        p->getX().operator=(utils::Vector<double, 3>(x_new));
                        Cells[i*ysize*zsize+zsize+l].addParticle(p);
                        Cells[i*ysize*zsize+(ysize - 1)*zsize+l].removeparticleAt(k);
                        k--;
                        break;
                    }
                    default: {
                        break;
                    }
                }
            }
        }
    }
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int i = 0; i < xsize; i++) {
        for (int j = 0; j < ysize; j++) {
            for (int k = 0; k < Cells[i*ysize*zsize+j*zsize].getElementNumber(); k++) {
                Particle* p = Cells[i*ysize*zsize+j*zsize].getParticleAt(k);
                switch (z1boundary) {
                    case 1: ///@brief remove particle if outflow
                    {
                        Cells[i*ysize*zsize+j*zsize].removeparticleAt(k);
                        k--;
                        break;
                    }
                    case 2: ///@brief if boundary reflective particle should not go out
                    {
                        break;
                    }
                    case 3: {///@brief set particle on the other side if periodic
                        double x_new[] = {p->getX().operator[](0), p->getX().operator[](1),
                                          (p->getX().operator[](2) + domain_size.operator[](2))};
                                           utils::Vector<int,3> a(0.0);
                        a[2]=1;
                        p->getPeriodicMovement().operator+(a);
                        p->getX().operator=(utils::Vector<double, 3>(x_new));
                        Cells[i*ysize*zsize+j*zsize+zsize - 2].addParticle(p);
                        Cells[i*ysize*zsize+j*zsize].removeparticleAt(k);
                        k--;
                        break;
                    }
                    default: {
                        break;
                    }
                }
            }
            for (int k = 0; k < Cells[i*ysize*zsize+j*zsize+zsize - 1].getElementNumber(); k++) {
                Particle *p = Cells[i*ysize*zsize+j*zsize+zsize - 1].getParticleAt(k);
                switch (z2boundary) {
                    case 1: ///@brief remove particle if outflow
                    {
                        Cells[i*ysize*zsize+j*zsize+zsize - 1].removeparticleAt(k);
                        k--;
                        break;
                    }
                    case 2: ///@brief if boundary reflective particle should not go out
                    {
                        break;
                    }
                    case 3:
                        ///@brief set particle on the other side if periodic
                    {
                        double x_new[] = {p->getX().operator[](0), p->getX().operator[](1),
                                          (p->getX().operator[](2) - domain_size.operator[](2))};
                                          utils::Vector<int,3> a(0.0);
                        a[2]=-1;
                        p->getPeriodicMovement().operator+(a);
                        p->getX().operator=(utils::Vector<double, 3>(x_new));
                        Cells[i*ysize*zsize+j*zsize+1].addParticle(p);
                        Cells[i*ysize*zsize+j*zsize+zsize - 1].removeparticleAt(k);
                        k--;
                        break;
                    }
                    default: {
                        break;
                    }
                }
            }
        }
    }
}

