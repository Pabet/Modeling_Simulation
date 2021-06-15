#include "particleGenerator2.h"
#include <list>
#include <cstdlib>

#include <vector>

using namespace std;
particleGenerator2::particleGenerator2(bool vo,bool vp,bool vrf,bool ho,bool hp,bool hrf, utils::Vector <double,3> domain, double r,bool mbd1) {
radius=r;
voutflow=vo;
vperiodic=vp;
vreflecting=vrf;
houtflow=ho;
hperiodic=hp;
hreflecting=hrf;
x=domain[0];
y=domain[1];
z=domain [2];
mbd=mbd1;
}

struct molecule {
    utils::Vector <double,3> force1;
    Particle c=Particle(1);
Particle *p=&c;

molecule(utils::Vector <double,3> f, Particle *d)
{
    force1.operator=(f);
    p=d;
}


};



void calculateForce2 (vector <molecule> &a, vector <molecule> &b,double sigma, double epsilon,bool mbd,double radius ){
    double l[3]={0,0,0};
    for (int i=0; i<a.size();i++) {

        Particle &p=*a[i].p;

        utils::Vector <double,3> c(l);
        for (int n=0; n<b.size();n++) {
                molecule &m = b[n];

                Particle &p2 = *b[n].p;
                if (mbd == true) {
                    if ((p.getX().operator-(p2.getX())).L2Norm()<=radius){
                    float k = (p.getX() - p2.getX()).L2Norm();
                    float l3 = sqrt(p.getsig()*p2.getsig())/ k;
                    float f1 = ((24 *(( p.geteps()+p2.geteps())/2)) / pow(k, 2)) * (pow(l3, 6) - 2 * pow(l3, 12));
                    utils::Vector<double, 3> b1 = (p2.getX().operator-(p.getX())).operator*(f1);
                    c.operator=(c.operator+(b1));
                    b[n].force1.operator=(b[n].force1.operator+(b1.operator*(-1)));
                }} else {
                    double scalar = (p.getM() * p2.getM()) / pow((p.getX().operator-(p2.getX())).L2Norm(), 3);
                    utils::Vector<double, 3> force = (p2.getX().operator-(p.getX())).operator*(scalar);
                    c.operator=(c.operator+(force));
                    b[n].force1.operator=(b[n].force1.operator+(force.operator*(-1)));
                }
                // b[n].force1.operator=(b[n].force1.operator+(b1.operator*(-1)));

            }

        a[i].force1.operator=(a[i].force1.operator+(c));

    }
}



void linkedCell (vector <vector<molecule>> &particles,double sigma,double epsilon, bool houtflow,bool hperiodic, bool hreflecting, bool voutflow,bool vperiodic, bool vreflecting, double a1, double b,double c, double radius,bool mbd,bool gravity,double g){
   int x= static_cast<int>( (a1 / radius));

   int y= static_cast<int>( (b / radius));

   int z= static_cast<int>( (c / radius));
    double l[3]={0,0,0};
    double factor=pow(2,1/6)*sigma;

    for (long unsigned int i=0;i<particles.size();i++) {

        for (long unsigned int k=0;k<particles[i].size();k++){

            Particle &p1=*particles[i][k].p;
            utils::Vector <double,3> c(l);

           for (long unsigned int p=k+1; p<particles[i].size();p++){
            
           Particle &p2=*particles[i][p].p;
           if (mbd==true){
               float k = (p1.getX() - p2.getX()).L2Norm();
               float l3=sqrt(p1.getsig()*p2.getsig())/ k;
               float f1 = (((24 * ( p1.geteps()+p2.geteps())/2)) / pow(k, 2)) * (pow(l3, 6) - 2 * pow(l3, 12));        ///@brief sums inside the cell over the vector of molecules
               utils::Vector<double, 3> b1 = (p2.getX().operator-(p1.getX())).operator*(f1);
               c.operator=(c.operator+(b1));
                particles[i][p].force1.operator=(particles[i][p].force1.operator+(b1.operator*(-1)));
           }
               else {double scalar = (p1.getM() * p2.getM()) / pow((p1.getX().operator-( p2.getX())).L2Norm(), 3);
            utils::Vector<double,3> force = (p2.getX().operator-(p1.getX())).operator*(scalar);
            c.operator=(c.operator+(force));
                   particles [i][p].force1.operator=(particles[i][p].force1.operator+(force.operator*(-1)));
               }
              // particles[i][p].force1.operator=(particles[i][p].force1.operator+(b1.operator*(-1)));

           }

          
            particles[i][k].force1.operator=(particles[i][k].force1.operator+(c));
        }

       if (!particles[i].empty()){

           Particle s=*particles[i][0].p;
            int position_x= static_cast<int>(floor ((s.getX()[0]) / radius));
int position_y= static_cast<int>( floor(s.getX()[1] / radius));
int position_z= static_cast<int>(floor (s.getX()[2] / radius));
            if ((position_x+1)<=(x-1) && (position_x+1)>=0){

                calculateForce2(particles[i],particles [position_x+1+position_y*x+position_z*y*x],sigma,epsilon,mbd,radius);
                if (z!=0 && position_z+1<z && position_z+1>=0){
                    calculateForce2(particles[i],particles [position_x+1+position_y*x+(position_z+1)*y*x],sigma,epsilon,mbd,radius);       ///@breif sums over the neighbouring cells
                }
            }
               if (vperiodic||hperiodic){
                if ((position_y!=0&& position_x!=0) || (position_y!=0 && position_x!=x)||( position_y!=y && position_x!=0) || (position_y!=y&& position_x!=x))
                {if (position_y==0 &&hperiodic){
                    calculateForce2(particles[i],particles [x*(y-1)+position_x-1+position_z*x*y],sigma,epsilon,mbd,radius);
                calculateForce2(particles[i],particles[x*(y-1)+position_x+position_z*x*y],sigma,epsilon, mbd,radius);
                calculateForce2(particles[i],particles[x*(y-1)+position_x-2+ position_z*x*y],sigma,epsilon,mbd,radius);

                }
                if (position_x==0&&vperiodic){
                    calculateForce2(particles[i],particles[(position_y+1)*x+x-1+position_z*x*y],sigma,epsilon,mbd,radius);
                    calculateForce2(particles [i],particles [i+x-1],sigma,epsilon,mbd,radius);
                // calculateForce2(particles[i],particles[(position_y-1)*x+x-1*position_z*x*y],sigma,epsilon,mbd,radius);

                }
                if (position_x==(x-1)&&vperiodic){
                    calculateForce2(particles[i],particles[position_y*x+position_z*x*y],sigma,epsilon,mbd,radius);
                }

                }
                if (position_x==0&&position_y==0&&vperiodic&&hperiodic){
                    calculateForce2(particles[i],particles[position_x+x-1+position_y*x+position_z*x*y],sigma,epsilon,mbd,radius);
                    calculateForce2(particles[i],particles[(position_y+1)*x+x-1+position_z*x*y],sigma,epsilon,mbd,radius);
                    calculateForce2(particles[i],particles[(y-1)*x+position_z*x*y],sigma,epsilon,mbd,radius);
                    calculateForce2(particles[i],particles[(y-1)*x+x-1+position_z*x*y],sigma,epsilon,mbd,radius);
                    calculateForce2(particles[i],particles[(y-1)*x+1+position_z*x*y],sigma,epsilon,mbd,radius);
                }
                if (position_x==(x-1)&&position_y==0&&vperiodic&&hperiodic){
                    calculateForce2(particles[i],particles[(y-1)*x+x-1+position_z*x*y],sigma,epsilon,mbd,radius);
                    calculateForce2(particles[i],particles[(y-1)*x+position_z*x*y],sigma,epsilon,mbd,radius);
                    calculateForce2(particles[i],particles[(y-1)*x+1+position_z*x*y],sigma,epsilon,mbd,radius);
                    calculateForce2(particles[i],particles[i+1],sigma,epsilon,mbd,radius);
                }
                if (position_y==(y-1)&&position_x==0&&vperiodic&&hperiodic){
                    calculateForce2(particles[i],particles[i+x-1],sigma,epsilon,mbd,radius);
                    //calculateForce2(particles[i],particles[(y-2)*x+x-1+position_z*x*y],sigma,epsilon,mbd,radius);

                }
                
                }
            
            
            if ((position_y+1)<=(y-1) && (position_y+1)>=0){
                calculateForce2(particles[i],particles [position_x+(position_y+1)*x+position_z*y*x],sigma,epsilon,mbd,radius);
                if (z!=0 && position_z+1<z && position_z+1>=0) {
                    calculateForce2(particles[i],particles [position_x+(position_y+1)*x+(position_z+1)*y*x],sigma,epsilon,mbd,radius);
                }
            }
            if ((position_y+1)<=(y-1) && (position_y+1)>=0 && (position_x-1)<=(x-1) && (position_x-1)>=0) {

                calculateForce2(particles[i],particles [position_x-1+(position_y+1)*x+position_z*y*x],sigma,epsilon,mbd,radius);
                if (z!=0 && position_z+1<z && position_z+1>=0) {
                    calculateForce2(particles[i],particles [position_x-1+(position_y+1)*x+(position_z+1)*y*x],sigma,epsilon,mbd,radius);
                }
            }
            if ((position_y+1)<=(y-1) && (position_y+1)>=0 && (position_x+1)<=(x-1) && (position_x+1)>=0) {

                calculateForce2(particles[i],particles [position_x+1+(position_y+1)*x+position_z*y*x],sigma,epsilon,mbd,radius);
                if (z!=0 && position_z+1<=(z-1) && position_z+1>=0) {
                    calculateForce2(particles[i],particles [position_x+1+(position_y+1)*x+(position_z+1)*y*x],sigma,epsilon,mbd,radius);
                }
            }

            if (z!=0) {
                if (position_z+1<z && position_z+1>=0) {
                    calculateForce2(particles[i],particles [position_x+position_y*x+(position_z+1)*y*x],sigma,epsilon,mbd,radius);
                    if (position_x-1>=0 && position_x-1<x){
                        calculateForce2(particles[i],particles [position_x-1+position_y*x+(position_z+1)*y*x],sigma,epsilon,mbd,radius);
                    }
                    if (position_y-1>=0 && position_y-1<y){
                        calculateForce2(particles[i],particles [position_x+(position_y-1)*x+(position_z+1)*y*x],sigma,epsilon,mbd,radius);     //ismeneno
                        if (position_x-1>=0 && position_x-1<x){
                            calculateForce2(particles[i],particles [position_x-1+(position_y-1)*x+(position_z+1)*y*x],sigma,epsilon,mbd,radius);
                        }
                        if (position_x+1<x && position_x+1 >=0){
                            calculateForce2(particles[i],particles [position_x+1+(position_y-1)*x+(position_z+1)*y*x],sigma,epsilon,mbd,radius);
                        }
                    }
                }

            }


        for (int inner_size=0; inner_size< particles[i].size();inner_size++) {
            Particle &p=*particles[i][inner_size].p; //!!
           if (vreflecting==true||hreflecting==true){         ///@brief which yy is the closest to the particle, adds repulsive force
                if (p.getX()[0]< factor&&vreflecting)
                {
                    double s[]={0,p.getX()[1],p.getX()[2]};
                    utils::Vector <double,3> p2=utils::Vector<double ,3>(s);
                    float k = (p.getX().operator-(p2)).L2Norm();
                    float l3 = sigma / k;
                    float f1 = ((24 * epsilon) / pow(k, 2)) * (pow(l3, 6) - 2 * pow(l3, 12));
                    particles[i][inner_size].force1.operator=(particles[i][inner_size].force1.operator+( (p2.operator-(p.getX())).operator*(f1)));}
                if ((a1-p.getX()[0])<factor&&vreflecting){
                    double s[]={a1,p.getX()[1],p.getX()[2]};
                    utils::Vector <double,3> p2=utils::Vector<double ,3>(s);
                    float k = (p.getX().operator-(p2)).L2Norm();
                    float l3 = sigma / k;
                    float f1 = ((24 * epsilon) / pow(k, 2)) * (pow(l3, 6) - 2 * pow(l3, 12));     //zagnat' v fct
                    particles[i][inner_size].force1.operator=(particles[i][inner_size].force1.operator+ ((p2.operator-(p.getX())).operator*(f1)));
                }
                if (p.getX()[1]<factor&&hreflecting) {
                    double s[]={p.getX()[0],0,p.getX()[2]};
                    utils::Vector <double,3> p2=utils::Vector<double ,3>(s);
                    float k = ((p.getX().operator-(p2)).L2Norm());
                    float l3 = (sigma / k);
                    float f1 = (((24 * epsilon) / pow(k, 2)) * (pow(l3, 6) - 2 * pow(l3, 12)));
                    particles[i][inner_size].force1.operator=(particles[i][inner_size].force1.operator+ ((p2.operator-(p.getX())).operator*(f1)));
                }
                if ((b-p.getX()[1])<factor&&hreflecting) {
                    double s[]={p.getX()[0],b,p.getX()[2]};
                    utils::Vector <double,3> p2=utils::Vector<double ,3>(s);
                    float k = (p.getX().operator-(p2)).L2Norm();
                    float l3 = sigma / k;
                    float f1 = ((24 * epsilon) / pow(k, 2)) * (pow(l3, 6) - 2 * pow(l3, 12));
                    particles[i][inner_size].force1.operator=(particles[i][inner_size].force1.operator+ ((p2.operator-(p.getX())).operator*(f1)));
                }

                if (p.getX()[2]<factor && c!=0&&vreflecting&&hreflecting) {
                    double s[]={p.getX()[0],p.getX()[1],0};
                    utils::Vector <double,3> p2=utils::Vector<double ,3>(s);
                    float k = (p.getX().operator-(p2)).L2Norm();
                    float l3 = sigma / k;
                    float f1 = ((24 * epsilon) / pow(k, 2)) * (pow(l3, 6) - 2 * pow(l3, 12));
                    particles[i][inner_size].force1.operator=(particles[i][inner_size].force1.operator+ ((p2.operator-(p.getX())).operator*(f1)));
                }
                if ((c-p.getX()[2])<factor && c!=0&&vreflecting&&hreflecting) {
                    double s[]={p.getX()[0],p.getX()[1],c};
                    utils::Vector <double,3> p2=utils::Vector<double ,3>(s);
                    float k = (p.getX().operator-(p2)).L2Norm();
                    float l3 = sigma / k;
                    float f1 = ((24 * epsilon) / pow(k, 2)) * (pow(l3, 6) - 2 * pow(l3, 12));
                    particles[i][inner_size].force1.operator=(particles[i][inner_size].force1.operator+ ((p2.operator-(p.getX())).operator*(f1)));
                }
            }
            p.getOldF().operator=(p.getF());
            p.getF().operator=(particles[i][inner_size].force1);
            if (gravity){
                p.getF()[1]=p.getF()[1]+p.getM()*g;
            }
        }



        }

    }





}







int &particleGenerator2::createMesh (std::list <Particle> &particles, double sigma, double epsilon,bool gravity, double g){
 std::list<Particle>::iterator iterator=particles.begin();
 int  size_x= static_cast<int>( (x / radius));
 int  size_y= static_cast<int>( (y / radius));
 int  size_z= static_cast<int> ((z / radius));
 double l[3]={0,0,0};
 utils::Vector <double,3> force2(l);
vector <vector <molecule>> a;

int size;
if (z!=0) size=(size_x) * (size_y) * (size_z);
else size=(size_x)*size_y;


vector <vector<molecule>> b(static_cast<unsigned long>(size));

while (iterator!=particles.end()){
 Particle &p=*iterator;
 if(voutflow||houtflow) {
     if ((((p.getX()[0] > x || p.getX()[0] < 0)&&voutflow) || ((p.getX()[1] > y || p.getX()[1] < 0)&&houtflow) ||
          ((p.getX()[2] > z || p.getX()[2] < 0) && z != 0&&voutflow&&houtflow))) {
         iterator = particles.erase(iterator);///@brief detect and erase halo particles
     }
 }else if(vperiodic||hperiodic){ ///@brief Move particles to other side
    while(p.getX()[0]>x&&vperiodic){
        double s[]={p.getX()[0]-x,p.getX()[1],p.getX()[2]};
        utils::Vector<double,3> temp=utils::Vector<double ,3>(s);
        p.getX().operator=(temp);
    }
    while(p.getX()[0]<0&&vperiodic){
        double s[]={p.getX()[0]+x,p.getX()[1],p.getX()[2]};
        utils::Vector<double,3> temp=utils::Vector<double ,3>(s);
        p.getX().operator=(temp);
    }
    while(p.getX()[1]>y&&hperiodic){
        double s[]={p.getX()[0],p.getX()[1]-y,p.getX()[2]};
        utils::Vector<double,3> temp=utils::Vector<double ,3>(s);
        p.getX().operator=(temp);
    }

    while(p.getX()[1]<0&&hperiodic){
        double s[]={p.getX()[0],p.getX()[1]+y,p.getX()[2]};
        utils::Vector<double,3> temp=utils::Vector<double ,3>(s);
        p.getX().operator=(temp);
    }
    if(z!=0){
        while(p.getX()[2]<0&&vperiodic&&hperiodic){
        double s[]={p.getX()[0],p.getX()[1],p.getX()[2]+z};
        utils::Vector<double,3> temp=utils::Vector<double ,3>(s);
        p.getX().operator=(temp);
        }
        while(p.getX()[2]>z&&vperiodic&&hperiodic){
            double s[]={p.getX()[0],p.getX()[1],p.getX()[2]-z};
            utils::Vector<double,3> temp=utils::Vector<double ,3>(s);
            p.getX().operator=(temp);
        }
    }
    ++iterator;
}else{Particle *pointer2=&(*iterator);
    if(pointer2!=NULL) {
        molecule d = {force2, pointer2};

        int position = static_cast<int>( floor(p.getX()[0] / radius) + (floor(p.getX()[1] / radius)) * size_x +
                                         (floor(p.getX()[2] / radius)) * size_x * size_y);

        if (position >= 0 && position < size)
        { b[position].push_back(d);  }                 ///@brief pack molecules into cells
    }
++iterator;
}


}
if(!b[0].empty())
linkedCell(b,sigma,epsilon,houtflow,hperiodic,hreflecting,voutflow,vperiodic,vreflecting,x, y, z,radius,mbd,gravity,g);     ///@brief calculate new force


}








