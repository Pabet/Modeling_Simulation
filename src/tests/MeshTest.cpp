//
// Created by Victor Stroescu on 12/01/2019.
//

#include "MeshTest.h"
#include "gtest/gtest.h"
#include "Particle.h"
#include "utils/Vector.h"
#include "Mesh.h"
#include "MaxwellBoltzmannDistribution.h"


TEST(initiateMesh,noparticles) {
    std::list<Particle> empty;
    double x[] = {1, 1, 0};
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1, 0, 0, 0, 0, 0, 0, 0);
    ASSERT_TRUE(c.particles.empty());
    ASSERT_TRUE(c.xsize == 3);
    ASSERT_TRUE(c.ysize == 3);
    ASSERT_TRUE(c.zsize == 3);
}TEST(initiateMesh,emptydomain)
{
    std::list<Particle> empty;
    double x[]={0,0,0};
    utils::Vector<double,3> d=utils::Vector<double ,3>(x);
    Mesh c=Mesh(false,0,1,d,empty,1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.empty());
    ASSERT_TRUE(c.xsize==2);
    ASSERT_TRUE(c.ysize==2);
    ASSERT_TRUE(c.zsize==3);
}
TEST(reasignMesh, eightParticlesinonecell)
{
    std::list<Particle> empty;
    double x[]={1,1,0};
    for(int i=1;i<=2;i++)
    {
        for(int j=1;j<=2;j++) {
            double c=double(i)/5;
            double d=double(j)/5;
            double cn=double(-i)/5;
            double dn=double(-j)/5;
            double position[] = {c,d,0};
            double positionn[] = {cn,dn,0};
            double v[]={0,0,0};
            Particle a = Particle(utils::Vector<double ,3>(position),utils::Vector<double ,3>(v), 1, 0, 0, 0 );
            Particle b = Particle(utils::Vector<double ,3>(positionn),utils::Vector<double ,3>(v), 1, 0, 0, 0 );
            empty.push_back(a);
            empty.push_back(b);
        }
    }
    utils::Vector<double,3> d=utils::Vector<double ,3>(x);
    Mesh c=Mesh(false,0,1,d,empty,1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size()==8);
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    ASSERT_TRUE(c.xsize==3);
    ASSERT_TRUE(c.ysize==3);
    ASSERT_TRUE(c.zsize==3);
    c.reassignMesh();
    EXPECT_TRUE(c.Cells[z1*y1+z1+1].getElementNumber()==8);
}
TEST(addandremove, addandremove32particles)
{
    particleContainer s=particleContainer();
    for(int i=0;i<5;i++) {
        for (int j =0; j <5; j++) {
            double c = double(i) / 5;
            double d = double(j) / 5;
            double cn = double(-i) / 5;
            double dn = double(-j) / 5;
            double positionpp[] = {c, d, 0};
            double positionnn[] = {cn, dn, 0};
            double positionpn[] = {c, dn, 0};
            double positionnp[] = {cn, d, 0};
            double v[] = {0, 0, 0};
            Particle* a =new Particle(utils::Vector<double, 3>(positionpp), utils::Vector<double, 3>(v), 1, 0, 0, 0);
            Particle an = Particle(utils::Vector<double, 3>(positionnn), utils::Vector<double, 3>(v), 1, 0, 0, 0);
            Particle b = Particle(utils::Vector<double, 3>(positionpn), utils::Vector<double, 3>(v), 1, 0, 0, 0);
            Particle bn = Particle(utils::Vector<double, 3>(positionnp), utils::Vector<double, 3>(v), 1, 0, 0, 0);
            s.addParticle(a);
            EXPECT_EQ(s.getElementNumber(),i*5+j+1);

        }
    }
    EXPECT_EQ(s.getElementNumber(),25);
    for(int i=24;i>=0;i--)
    {
        s.removeparticleAt(i);
    }
    for(int i=0;i<6;i++) {
        for (int j =0; j <6; j++) {
            double c = double(i) / 5;
            double d = double(j) / 5;
            double cn = double(-i) / 5;
            double dn = double(-j) / 5;
            double positionpp[] = {c, d, 0};
            double positionnn[] = {cn, dn, 0};
            double positionpn[] = {c, dn, 0};
            double positionnp[] = {cn, d, 0};
            double v[] = {0, 0, 0};
            Particle* a =new Particle(utils::Vector<double, 3>(positionpp), utils::Vector<double, 3>(v), 1, 0, 0, 0);
            Particle an = Particle(utils::Vector<double, 3>(positionnn), utils::Vector<double, 3>(v), 1, 0, 0, 0);
            Particle b = Particle(utils::Vector<double, 3>(positionpn), utils::Vector<double, 3>(v), 1, 0, 0, 0);
            Particle bn = Particle(utils::Vector<double, 3>(positionnp), utils::Vector<double, 3>(v), 1, 0, 0, 0);
            s.addParticle(a);
            EXPECT_EQ(s.getElementNumber(),i*6+j+1);

        }
    }
    for(int i=35;i>=0;i--)
    {
        s.removeparticleAt(i);
    }
}
TEST(reasignMesh, sixteenParticlesin4cells)
{
    std::list<Particle> empty;
    double x[]={1.0,1.0,0.0};
    for(int i=1;i<=2;i++)
    {
        for(int j=1;j<=2;j++) {
            double c=double(i)/5;
            double d=double(j)/5;
            double cn=double(-i)/5;
            double dn=double(-j)/5;
            double positionpp[] = {c,d,0};
            double positionnn[] = {cn,dn,0};
            double positionpn[]={c,dn,0};
            double positionnp[]={cn,d,0};
            double v[]={0,0,0};
            Particle a = Particle(utils::Vector<double ,3>(positionpp),utils::Vector<double ,3>(v), 1, 0, 0, 0 );
            Particle an = Particle(utils::Vector<double ,3>(positionnn),utils::Vector<double ,3>(v), 1, 0, 0, 0 );
            Particle b = Particle(utils::Vector<double ,3>(positionpn),utils::Vector<double ,3>(v), 1, 0, 0, 0 );
            Particle bn = Particle(utils::Vector<double ,3>(positionnp),utils::Vector<double ,3>(v), 1, 0, 0, 0 );
            empty.push_back(a);
            empty.push_back(b);
            empty.push_back(an);
            empty.push_back(bn);
        }
    }
    utils::Vector<double,3> d=utils::Vector<double ,3>(x);
    Mesh c=Mesh(false,0,0.5,d,empty,1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size()==16);
    ASSERT_TRUE(c.xsize==4);
    ASSERT_TRUE(c.ysize==4);
    ASSERT_TRUE(c.zsize==3);
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    c.reassignMesh();
    EXPECT_TRUE(c.Cells[z1*y1+z1+1].getElementNumber()==4);
    EXPECT_TRUE(c.Cells[2*z1*y1+2*z1+1].getElementNumber()==4);
    EXPECT_TRUE(c.Cells[z1*y1+2*z1+1].getElementNumber()==4);
    EXPECT_TRUE(c.Cells[2*z1*y1+z1+1].getElementNumber()==4);
}
TEST(reasignMesh, nineParticlesinonecell)
{
    std::list<Particle> empty;
    double x[]={1,1,0};
    for(int i=-1;i<=1;i++)
    {
        for(int j=-1;j<=1;j++) {
            double c=double(i)/3;
            double d=double(j)/3;
            double position[] = {c,d,0};
            double v[]={0,0,0};
            Particle a = Particle(utils::Vector<double ,3>(position),utils::Vector<double ,3>(v), 1, 0, 0, 0 );
            empty.push_back(a);
        }
    }
    utils::Vector<double,3> d=utils::Vector<double ,3>(x);
    Mesh c=Mesh(false,0,1,d,empty,1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size()==9);
    ASSERT_TRUE(c.xsize==3);
    ASSERT_TRUE(c.ysize==3);
    ASSERT_TRUE(c.zsize==3);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    EXPECT_TRUE(c.Cells[z1*y1+z1+1].getElementNumber()==9);
}
TEST(reasignMesh, nineParticlesinninecells)
{
    std::list<Particle> empty;
    double x[]={3,3,0};
    for(int i=-1;i<=1;i++)
    {
        for(int j=-1;j<=1;j++) {
            double c=i;
            double d=j;
            double position[] = {c,d,0};
            double v[]={0,0,0};
            Particle a = Particle(utils::Vector<double ,3>(position),utils::Vector<double ,3>(v), 1, 0, 0, 0 );
            empty.push_back(a);
        }
    }
    utils::Vector<double,3> d=utils::Vector<double ,3>(x);
    Mesh c=Mesh(false,0,1,d,empty,1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size()==9);
    ASSERT_TRUE(c.xsize==5);
    ASSERT_TRUE(c.ysize==5);
    ASSERT_TRUE(c.zsize==3);
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    c.reassignMesh();
    for(int i=1;i<=3;i++)
        for(int j=1;j<=3;j++)
            EXPECT_TRUE(c.Cells[i*z1*y1+j*z1+1].getElementNumber()==1);
}
TEST(reasignMesh, twentysevenParticlesinonecell3D)
{
    std::list<Particle> empty;
    double x[]={1,1,1};
    for(int i=-1;i<=1;i++)
    {
        for(int j=-1;j<=1;j++) {
            for (int l = -1; l <= 1; l++) {
                double c = (double)i / 3;
                double d = (double)j / 3;
                double e = (double)l / 3;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double ,3>(position),utils::Vector<double ,3>(v), 1, 0, 0, 0 );
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double,3> d=utils::Vector<double ,3>(x);
    Mesh c=Mesh(false,0,1,d,empty,1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size()==27);
    ASSERT_TRUE(c.xsize==3);
    ASSERT_TRUE(c.ysize==3);
    ASSERT_TRUE(c.zsize==3);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    EXPECT_TRUE(c.Cells[z1*y1+z1+1].getElementNumber()==27);
}
TEST(reasignMesh, onehundredtwentyoneParticlesinonecell)
{
    std::list<Particle> empty;
    double x[]={1,1,0};
    for(int i=-5;i<=5;i++)
    {
        for(int j=-5;j<=5;j++) {
            double c=(double)i/11;
            double d=(double)j/11;
            double position[] = {c,d,0};
            double v[]={0,0,0};
            Particle a = Particle(utils::Vector<double ,3>(position),utils::Vector<double ,3>(v), 1, 0, 0, 0 );
            empty.push_back(a);
        }
    }
    utils::Vector<double,3> d=utils::Vector<double ,3>(x);
    Mesh c=Mesh(false,0,1,d,empty,1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size()==121);
    ASSERT_TRUE(c.xsize==3);
    ASSERT_TRUE(c.ysize==3);
    ASSERT_TRUE(c.zsize==3);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    EXPECT_TRUE(c.Cells[z1*y1+z1+1].getElementNumber()==121);
}

TEST(reasignMesh, onehundredtwentyoneParticlesinonehundredtwentyonecells)
{
    std::list<Particle> empty;
    double x[]={11,11,0};
    for(int i=-5;i<=5;i++)
    {
        for(int j=-5;j<=5;j++) {

                double c = i;
                double d = j;
                double position[] = {c, d, 0};
                double v[] = {0, 0, 0};
            Particle a = Particle(utils::Vector<double ,3>(position),utils::Vector<double ,3>(v), 1, 0, 0, 0 );
                empty.push_back(a);

        }
    }
    utils::Vector<double,3> d=utils::Vector<double ,3>(x);
    Mesh c=Mesh(false,0,1,d,empty,1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size()==121);
    ASSERT_TRUE(c.xsize==13);
    ASSERT_TRUE(c.ysize==13);
    ASSERT_TRUE(c.zsize==3);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    for(int i=1;i<=11;i++) {
        for (int j = 1; j <= 11; j++) {
            EXPECT_TRUE(c.Cells[i*z1*y1+j*z1+1].getElementNumber() == 1);
        }
    }
    for(int i=1;i<=11;i++) {
            EXPECT_FALSE(c.Cells[i*z1*y1+1].getElementNumber() != 0);
            EXPECT_FALSE(c.Cells[z1*i+1].getElementNumber() != 0);
            EXPECT_FALSE(c.Cells[12*y1*z1+z1*i+1].getElementNumber() != 0);
            EXPECT_FALSE(c.Cells[i*y1*z1+z1*12+1].getElementNumber() != 0);
    }
}
TEST(reasignMesh, onethreethreeoneParticlesinonethreethreeonecells3d)
{
    std::list<Particle> empty;
    double x[]={11,11,11};
    for(int i=-5;i<=5;i++)
    {
        for(int j=-5;j<=5;j++) {
            for (int l = -5; l <= 5; l++) {
                double c = i;
                double d = j;
                double e = l;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double ,3>(position),utils::Vector<double ,3>(v), 1, 0, 0, 0 );
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double,3> d=utils::Vector<double ,3>(x);
    Mesh c=Mesh(false,0,1,d,empty,1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size()==1331);
    ASSERT_TRUE(c.xsize==13);
    ASSERT_TRUE(c.ysize==13);
    ASSERT_TRUE(c.zsize==13);
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    c.reassignMesh();
    for(int i=1;i<=11;i++) {
        for (int j = 1; j <= 11; j++) {
            for(int l=1;l<=11;l++) {
                EXPECT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    for(int i=1;i<=11;i++) {
        for(int j=1;j<11;j++) {
            EXPECT_FALSE(c.Cells[i*z1*y1+1].getElementNumber() != 0);
            EXPECT_FALSE(c.Cells[z1*i+1].getElementNumber() != 0);
            EXPECT_FALSE(c.Cells[12*y1*z1+z1*i+1].getElementNumber() != 0);
            EXPECT_FALSE(c.Cells[i*y1*z1+z1*12+1].getElementNumber() != 0);
            EXPECT_FALSE(c.Cells[i*y1*z1+z1*j].getElementNumber() != 0);
            EXPECT_FALSE(c.Cells[i*y1*z1+z1*j+12].getElementNumber() != 0);
        }

    }
}
TEST(testforce,testforce1particlealone)
{
    std::list<Particle> empty;
    double x[]={1,1,0};
    utils::Vector<double,3> d=utils::Vector<double ,3>(x);
    double position[] = {0, 0, 0};
    double v[] = {0, 0, 0};
    Particle a = Particle(utils::Vector<double ,3>(position),utils::Vector<double ,3>(v), 1, 0, 0, 0 );
    empty.push_back(a);
    Mesh c=Mesh(false,0,1,d,empty,1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size()==1);
    ASSERT_TRUE(c.xsize==3);
    ASSERT_TRUE(c.ysize==3);
    ASSERT_TRUE(c.zsize==3);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    ASSERT_TRUE(c.Cells[y1*z1+z1+1].getElementNumber()==1);
    c.calculateF();
    ASSERT_TRUE(a.getF().operator==(utils::Vector<double ,3>(0.0)));
    ASSERT_TRUE(a.getOldF().operator==(utils::Vector<double ,3>(0.0)));
    empty.clear();
}
TEST(testforce,testforce2particlesinsamecellatdistancegreaterthan6thsquarerootofsigma)
{
    double sigma=1;
    double epsilon=5;
    double m=1;
    std::list<Particle> empty;
    double x[]={2,2,0};
    utils::Vector<double,3> d=utils::Vector<double ,3>(x);
    double positionp[] = {0.5, 0.5, 0};
    double positionn[] = {-0.5, -0.5, 0};
    double v[] = {0, 0, 0};
    Particle a = Particle(utils::Vector<double ,3>(positionp),utils::Vector<double ,3>(v), m, 0, epsilon, sigma );
    empty.push_front(a);
    Particle b = Particle(utils::Vector<double ,3>(positionn),utils::Vector<double ,3>(v), m, 0, epsilon, sigma );
    empty.push_front(b);
    Mesh c=Mesh(false,0,2.0,d,empty,1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size()==2);
    ASSERT_TRUE(c.xsize==3);
    ASSERT_TRUE(c.ysize==3);
    ASSERT_TRUE(c.zsize==3);
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    c.reassignMesh();
    ASSERT_TRUE(c.Cells[y1*z1+z1+1].getElementNumber()==2);
    ASSERT_TRUE((c.Cells[y1*z1+z1+1].getParticleAt(0)->getX().operator-(c.Cells[y1*z1+z1+1].getParticleAt(1)->getX())).L2Norm()>(pow(2,1/6)*sigma));
    c.calculateF();
    ASSERT_TRUE(c.Cells[y1*z1+z1+1].getParticleAt(0)->getF()==c.Cells[y1*z1+z1+1].getParticleAt(1)->getF().operator*(-1));
    EXPECT_TRUE(c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](0)>0&&c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](1)>0&&c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](2)==0);
    EXPECT_TRUE(c.Cells[y1*z1+z1+1].getParticleAt(1)->getF().operator[](0)<0&&c.Cells[y1*z1+z1+1].getParticleAt(1)->getF().operator[](1)<0&&c.Cells[y1*z1+z1+1].getParticleAt(1)->getF().operator[](2)==0);

}
TEST(testforce,testforce2particlesinsamecellatdistancelowerthan6thsquarerootofsigma)
{
    double sigma=1;
    double epsilon=5;
    double m=1;
    std::list<Particle> empty;
    double x[]={2,2,0};
    utils::Vector<double,3> d=utils::Vector<double ,3>(x);
    double positionp[] = {0.2, 0.2, 0};
    double positionn[] = {-0.2, -0.2, 0};
    double v[] = {0, 0, 0};
    Particle a = Particle(utils::Vector<double ,3>(positionp),utils::Vector<double ,3>(v), m, 0, epsilon, sigma );
    empty.push_front(a);
    Particle b = Particle(utils::Vector<double ,3>(positionn),utils::Vector<double ,3>(v), m, 0, epsilon, sigma );
    empty.push_front(b);
    Mesh c=Mesh(false,0,2.0,d,empty,1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size()==2);
    ASSERT_TRUE(c.xsize==3);
    ASSERT_TRUE(c.ysize==3);
    ASSERT_TRUE(c.zsize==3);
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    c.reassignMesh();
    ASSERT_TRUE(c.Cells[y1*z1+z1+1].getElementNumber()==2);
    ASSERT_TRUE((c.Cells[y1*z1+z1+1].getParticleAt(0)->getX().operator-(c.Cells[y1*z1+z1+1].getParticleAt(1)->getX())).L2Norm()<(pow(2,1/6)*sigma));
    c.calculateF();
    ASSERT_TRUE(c.Cells[y1*z1+z1+1].getParticleAt(0)->getF()==c.Cells[y1*z1+z1+1].getParticleAt(1)->getF().operator*(-1));
    EXPECT_TRUE(c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](0)<0&&c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](1)<0&&c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](2)==0);
    EXPECT_TRUE(c.Cells[y1*z1+z1+1].getParticleAt(1)->getF().operator[](0)>0&&c.Cells[y1*z1+z1+1].getParticleAt(1)->getF().operator[](1)>0&&c.Cells[y1*z1+z1+1].getParticleAt(1)->getF().operator[](2)==0);

}
TEST(testforce,testforce2particlesinsamecellatdistanceequaltothan6thsquarerootofsigma)
{
    double sigma=1;
    double epsilon=5;
    double m=1;
    std::list<Particle> empty;
    double x[]={2,2,0};
    utils::Vector<double,3> d=utils::Vector<double ,3>(x);
    double positionp[] = {((pow(2,2.0/3)*sigma)/4), ((pow(2,2.0/3)*sigma)/4), 0};
    double positionn[] = {-((pow(2,2.0/3)*sigma)/4), -((pow(2,2.0/3)*sigma)/4), 0};
    double v[] = {0, 0, 0};
    Particle a = Particle(utils::Vector<double ,3>(positionp),utils::Vector<double ,3>(v), m, 0, epsilon, sigma );
    empty.push_front(a);
    Particle b = Particle(utils::Vector<double ,3>(positionn),utils::Vector<double ,3>(v), m, 0, epsilon, sigma );
    empty.push_front(b);
    Mesh c=Mesh(false,0,2.0,d,empty,1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size()==2);
    ASSERT_TRUE(c.xsize==3);
    ASSERT_TRUE(c.ysize==3);
    ASSERT_TRUE(c.zsize==3);
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    c.reassignMesh();
    ASSERT_TRUE(c.Cells[y1*z1+z1+1].getElementNumber()==2);
    ASSERT_DOUBLE_EQ((c.Cells[y1*z1+z1+1].getParticleAt(0)->getX().operator-(c.Cells[y1*z1+z1+1].getParticleAt(1)->getX())).L2Norm(),(pow(2,1.0/6)*sigma));
    c.calculateF();
    ASSERT_TRUE(c.Cells[y1*z1+z1+1].getParticleAt(0)->getF()==c.Cells[y1*z1+z1+1].getParticleAt(1)->getF().operator*(-1));
    EXPECT_NEAR(c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](0),0,4*pow(10.0,-14));
    EXPECT_NEAR(c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](1),0,4*pow(10.0,-14));
    EXPECT_NEAR(c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](2),0,4*pow(10.0,-14));
    EXPECT_NEAR(c.Cells[y1*z1+z1+1].getParticleAt(1)->getF().operator[](0),0,4*pow(10.0,-14));
    EXPECT_NEAR(c.Cells[y1*z1+z1+1].getParticleAt(1)->getF().operator[](1),0,4*pow(10.0,-14));
    EXPECT_NEAR(c.Cells[y1*z1+z1+1].getParticleAt(1)->getF().operator[](2),0,4*pow(10.0,-14));
}TEST(testforce,testforce2particlesindifferentcellsatdistancegreaterthan6thsquarerootofsigma)
{
    double sigma=1;
    double epsilon=5;
    double m=1;
    std::list<Particle> empty;
    double x[]={3,3,0};
    utils::Vector<double,3> d=utils::Vector<double ,3>(x);
    double positionp[] = {0.5, 0.5, 0};
    double positionn[] = {-0.5, -0.5, 0};
    double v[] = {0, 0, 0};
    Particle a = Particle(utils::Vector<double ,3>(positionp),utils::Vector<double ,3>(v), m, 0, epsilon, sigma );
    empty.push_front(a);
    Particle b = Particle(utils::Vector<double ,3>(positionn),utils::Vector<double ,3>(v), m, 0, epsilon, sigma );
    empty.push_front(b);
    Mesh c=Mesh(false,0,1.5,d,empty,1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size()==2);
    ASSERT_TRUE(c.xsize==4);
    ASSERT_TRUE(c.ysize==4);
    ASSERT_TRUE(c.zsize==3);
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    c.reassignMesh();
    ASSERT_TRUE(c.Cells[y1*z1+z1+1].getElementNumber()==1);
    ASSERT_TRUE(c.Cells[2*y1*z1+2*z1+1].getElementNumber()==1);
    ASSERT_TRUE((c.Cells[y1*z1+z1+1].getParticleAt(0)->getX().operator-(c.Cells[2*y1*z1+2*z1+1].getParticleAt(0)->getX())).L2Norm()>(pow(2,1/6)*sigma));
    c.calculateF();
    ASSERT_TRUE(c.Cells[y1*z1+z1+1].getParticleAt(0)->getF()==c.Cells[2*y1*z1+2*z1+1].getParticleAt(0)->getF().operator*(-1));
    EXPECT_TRUE(c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](0)>0&&c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](1)>0&&c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](2)==0);
    EXPECT_TRUE(c.Cells[2*y1*z1+2*z1+1].getParticleAt(0)->getF().operator[](0)<0&&c.Cells[2*y1*z1+2*z1+1].getParticleAt(0)->getF().operator[](1)<0&&c.Cells[2*y1*z1+2*z1+1].getParticleAt(0)->getF().operator[](2)==0);

}
TEST(testforce,testforce2particlesindifferentcellsatdistancelowerthan6thsquarerootofsigma)
{
    double sigma=1;
    double epsilon=5;
    double m=1;
    std::list<Particle> empty;
    double x[]={2,2,0};
    utils::Vector<double,3> d=utils::Vector<double ,3>(x);
    double positionp[] = {0.2, 0.2, 0};
    double positionn[] = {-0.2, -0.2, 0};
    double v[] = {0, 0, 0};
    Particle a = Particle(utils::Vector<double ,3>(positionp),utils::Vector<double ,3>(v), m, 0, epsilon, sigma );
    empty.push_front(a);
    Particle b = Particle(utils::Vector<double ,3>(positionn),utils::Vector<double ,3>(v), m, 0, epsilon, sigma );
    empty.push_front(b);
    Mesh c=Mesh(false,0,1.0,d,empty,1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size()==2);
    ASSERT_TRUE(c.xsize==4);
    ASSERT_TRUE(c.ysize==4);
    ASSERT_TRUE(c.zsize==3);
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    c.reassignMesh();
    ASSERT_TRUE(c.Cells[y1*z1+z1+1].getElementNumber()==1);
    ASSERT_TRUE(c.Cells[2*y1*z1+2*z1+1].getElementNumber()==1);
    ASSERT_TRUE((c.Cells[y1*z1+z1+1].getParticleAt(0)->getX().operator-(c.Cells[2*y1*z1+2*z1+1].getParticleAt(0)->getX())).L2Norm()<(pow(2,1/6)*sigma));
    c.calculateF();
    ASSERT_TRUE(c.Cells[y1*z1+z1+1].getParticleAt(0)->getF()==c.Cells[2*y1*z1+2*z1+1].getParticleAt(0)->getF().operator*(-1));
    EXPECT_TRUE(c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](0)<0&&c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](1)<0&&c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](2)==0);
    EXPECT_TRUE(c.Cells[2*y1*z1+2*z1+1].getParticleAt(0)->getF().operator[](0)>0&&c.Cells[2*y1*z1+2*z1+1].getParticleAt(0)->getF().operator[](1)>0&&c.Cells[2*y1*z1+2*z1+1].getParticleAt(0)->getF().operator[](2)==0);

}
TEST(testforce,testforce2particlesindifferentcellsatdistanceequaltothan6thsquarerootofsigma)
{
    double sigma=1;
    double epsilon=5;
    double m=1;
    std::list<Particle> empty;
    double x[]={2,2,0};
    utils::Vector<double,3> d=utils::Vector<double ,3>(x);
    double positionp[] = {((pow(2,2.0/3)*sigma)/4), ((pow(2,2.0/3)*sigma)/4), 0};
    double positionn[] = {-((pow(2,2.0/3)*sigma)/4), -((pow(2,2.0/3)*sigma)/4), 0};
    double v[] = {0, 0, 0};
    Particle a = Particle(utils::Vector<double ,3>(positionp),utils::Vector<double ,3>(v), m, 0, epsilon, sigma );
    empty.push_front(a);
    Particle b = Particle(utils::Vector<double ,3>(positionn),utils::Vector<double ,3>(v), m, 0, epsilon, sigma );
    empty.push_front(b);
    Mesh c=Mesh(false,0,1.0,d,empty,1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size()==2);
    ASSERT_TRUE(c.xsize==4);
    ASSERT_TRUE(c.ysize==4);
    ASSERT_TRUE(c.zsize==3);
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    c.reassignMesh();
    ASSERT_TRUE(c.Cells[y1*z1+z1+1].getElementNumber()==1);
    ASSERT_TRUE(c.Cells[2*y1*z1+2*z1+1].getElementNumber()==1);
    ASSERT_DOUBLE_EQ((c.Cells[y1*z1+z1+1].getParticleAt(0)->getX().operator-(c.Cells[2*y1*z1+2*z1+1].getParticleAt(0)->getX())).L2Norm(),(pow(2,1.0/6)*sigma));
    c.calculateF();
    ASSERT_TRUE(c.Cells[y1*z1+z1+1].getParticleAt(0)->getF()==c.Cells[2*y1*z1+2*z1+1].getParticleAt(0)->getF().operator*(-1));
    EXPECT_NEAR(c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](0),0,4*pow(10.0,-14));
    EXPECT_NEAR(c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](1),0,4*pow(10.0,-14));
    EXPECT_NEAR(c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](2),0,4*pow(10.0,-14));
    EXPECT_NEAR(c.Cells[2*y1*z1+2*z1+1].getParticleAt(0)->getF().operator[](0),0,4*pow(10.0,-14));
    EXPECT_NEAR(c.Cells[2*y1*z1+2*z1+1].getParticleAt(0)->getF().operator[](1),0,4*pow(10.0,-14));
    EXPECT_NEAR(c.Cells[2*y1*z1+2*z1+1].getParticleAt(0)->getF().operator[](2),0,4*pow(10.0,-14));
}
TEST(forceCalculation, summofforcesofallpariclesequals0) {
    double sigma = 1;
    double epsilon = 5;
    double m = 1;
    std::list<Particle> empty;
    double x[] = {11, 11, 11};
    for (int i = -5; i <= 5; i++) {
        for (int j = -5; j <= 5; j++) {
            for (int l = -5; l <= 5; l++) {
                double c = i;
                double d = j;
                double e = l;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), m, 0, epsilon, sigma);
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size()==1331);
    ASSERT_TRUE(c.xsize==13);
    ASSERT_TRUE(c.ysize==13);
    ASSERT_TRUE(c.zsize==13);
    c.reassignMesh();
    c.calculateF();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    utils::Vector<double,3>summ=utils::Vector<double,3>(0.0);
    for(int i=1;i<=11;i++) {
        for (int j = 1; j <= 11; j++) {
            for(int l=1;l<=11;l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
                if(l==1)
                {
                    EXPECT_TRUE(c.Cells[i*y1*z1+j*z1+l].getParticleAt(0)->getF().operator[](2)==-120);
                }
                if(l==11)
                {
                    EXPECT_TRUE(c.Cells[i*y1*z1+j*z1+l].getParticleAt(0)->getF().operator[](2)==120);
                }
                if(l!=11&&l!=1)
                {
                    EXPECT_TRUE(c.Cells[i*y1*z1+j*z1+l].getParticleAt(0)->getF().operator[](2)==0);
                }
                summ.operator=(summ.operator+(c.Cells[i*y1*z1+j*z1+l].getParticleAt(0)->getF()));
           }
        }
    }
    EXPECT_TRUE(summ.L2Norm()==0);
}
TEST(hallocellsdecisions,x1outflow) {

    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    for (int i = -2; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            for (int l = -1; l <= 1; l++) {
                double c = i;
                double d = j;
                double e = l;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), 1, 0, 0, 0);
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size() == 36);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    for (int i = 0; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    c.x1boundary = 1;
    c.halloCellsDecisions();
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    for (int j=1;j<=3;j++){
        for(int l=1;l<=3;l++)
        {
            EXPECT_TRUE(c.Cells[j*z1+l].getElementNumber()==0);
        }
    }
}
TEST(hallocellsdecisions,x2outflow) {

    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    for (int i = -1; i <= 2; i++) {
        for (int j = -1; j <= 1; j++) {
            for (int l = -1; l <= 1; l++) {
                double c = i;
                double d = j;
                double e = l;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), 1, 0, 0, 0);
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size() == 36);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    for (int i = 1; i <= 4; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    c.x2boundary = 1;
    c.halloCellsDecisions();
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    for (int j=1;j<=3;j++){
        for(int l=1;l<=3;l++)
        {
            EXPECT_TRUE(c.Cells[4*y1*z1+j*z1+l].getElementNumber()==0);
        }
    }
}
TEST(hallocellsdecisions,y1outflow) {

    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    for (int i = -1; i <= 1; i++) {
        for (int j = -2; j <= 1; j++) {
            for (int l = -1; l <= 1; l++) {
                double c = i;
                double d = j;
                double e = l;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), 1, 0, 0, 0);
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size() == 36);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    for (int i = 1; i <= 3; i++) {
        for (int j = 0; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    c.y1boundary = 1;
    c.halloCellsDecisions();
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    for (int i=1;i<=3;i++){
        for(int l=1;l<=3;l++)
        {
            EXPECT_TRUE(c.Cells[i*y1*z1+l].getElementNumber()==0);
        }
    }
}
TEST(hallocellsdecisions,y2outflow) {

    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 2; j++) {
            for (int l = -1; l <= 1; l++) {
                double c = i;
                double d = j;
                double e = l;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), 1, 0, 0, 0);
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size() == 36);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 4; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    c.y2boundary = 1;
    c.halloCellsDecisions();
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    for (int i=1;i<=3;i++){
        for(int l=1;l<=3;l++)
        {
            EXPECT_TRUE(c.Cells[i*y1*z1+4*z1+l].getElementNumber()==0);
        }
    }
}
TEST(hallocellsdecisions,z1outflow) {

    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            for (int l = -2; l <= 1; l++) {
                double c = i;
                double d = j;
                double e = l;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), 1, 0, 0, 0);
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size() == 36);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 0; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    c.z1boundary = 1;
    c.halloCellsDecisions();
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    for (int i=1;i<=3;i++){
        for(int j=1;j<=3;j++)
        {
            EXPECT_TRUE(c.Cells[i*y1*z1+j*z1].getElementNumber()==0);
        }
    }
}
TEST(hallocellsdecisions,z2outflow) {

    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            for (int l = -1; l <= 2; l++) {
                double c = i;
                double d = j;
                double e = l;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), 1, 0, 0, 0);
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size() == 36);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 4; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    c.z2boundary = 1;
    c.halloCellsDecisions();
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    for (int i=1;i<=3;i++){
        for(int j=1;j<=3;j++)
        {
            EXPECT_TRUE(c.Cells[i*y1*z1+j*z1+4].getElementNumber()==0);
        }
    }
}
TEST(hallocellsdecisions,alloutflow) {

    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    for (int i = -2; i <= 2; i++) {
        for (int j = -2; j <= 2; j++) {
            for (int l = -2; l <= 2; l++) {
                double c = i;
                double d = j;
                double e = l;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), 1, 0, 0, 0);
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size() == 125);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    for (int i = 0; i <= 4; i++) {
        for (int j = 0; j <= 4; j++) {
            for (int l = 0; l <= 4; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    c.x1boundary = 1;
    c.x2boundary = 1;
    c.y1boundary = 1;
    c.y2boundary = 1;
    c.z1boundary = 1;
    c.z2boundary = 1;
    c.halloCellsDecisions();
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    for (int i=1;i<=3;i++){
        for(int j=1;j<=3;j++)
        {
            EXPECT_TRUE(c.Cells[i*y1*z1+j*z1+4].getElementNumber()==0);
            EXPECT_TRUE(c.Cells[i*y1*z1+j*z1].getElementNumber()==0);
            EXPECT_TRUE(c.Cells[i*y1*z1+j].getElementNumber()==0);
            EXPECT_TRUE(c.Cells[i*y1*z1+4*z1+j].getElementNumber()==0);
            EXPECT_TRUE(c.Cells[i*z1+j].getElementNumber()==0);
            EXPECT_TRUE(c.Cells[4*y1*z1+i*z1+j].getElementNumber()==0);
        }
    }
}
TEST(hallocellsdecisions,x1reflective) {
    double sigma = 1;
    double epsilon = 5;
    double m = 1;
    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    double position[] = {-1, 0, 0};
    double v[] = {0, 0, 0};
    Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), m, 0, epsilon, sigma);
    empty.push_back(a);
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    c.reassignMesh();
    ASSERT_TRUE(c.particles.size() == 1);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.x1boundary=2;
    c.calculateF();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    ASSERT_TRUE(c.particles.size() == 1);
    EXPECT_TRUE(c.Cells[y1*z1+2*z1+2].getElementNumber()==1);
    EXPECT_TRUE(c.Cells[y1*z1+2*z1+2].getParticleAt(0)->getF().operator[](0)==120);
    EXPECT_TRUE(c.Cells[y1*z1+2*z1+2].getParticleAt(0)->getF().operator[](1)==0);
    EXPECT_TRUE(c.Cells[y1*z1+2*z1+2].getParticleAt(0)->getF().operator[](2)==0);
}
TEST(hallocellsdecisions,x2reflective) {
    double sigma = 1;
    double epsilon = 5;
    double m = 1;
    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    double position[] = {1, 0, 0};
    double v[] = {0, 0, 0};
    Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), m, 0, epsilon, sigma);
    empty.push_back(a);
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    c.reassignMesh();
    ASSERT_TRUE(c.particles.size() == 1);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.x2boundary=2;
    c.calculateF();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    ASSERT_TRUE(c.particles.size() == 1);
    EXPECT_TRUE(c.Cells[3*y1*z1+2*z1+2].getElementNumber()==1);
    EXPECT_TRUE(c.Cells[3*y1*z1+2*z1+2].getParticleAt(0)->getF().operator[](0)==-120);
    EXPECT_TRUE(c.Cells[3*y1*z1+2*z1+2].getParticleAt(0)->getF().operator[](1)==0);
    EXPECT_TRUE(c.Cells[3*y1*z1+2*z1+2].getParticleAt(0)->getF().operator[](2)==0);
}
TEST(hallocellsdecisions,y1reflective) {
    double sigma = 1;
    double epsilon = 5;
    double m = 1;
    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    double position[] = {0, -1, 0};
    double v[] = {0, 0, 0};
    Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), m, 0, epsilon, sigma);
    empty.push_back(a);
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    c.reassignMesh();
    ASSERT_TRUE(c.particles.size() == 1);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.y1boundary=2;
    c.calculateF();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    ASSERT_TRUE(c.particles.size() == 1);
    EXPECT_TRUE(c.Cells[2*y1*z1+z1+2].getElementNumber()==1);
    EXPECT_TRUE(c.Cells[2*y1*z1+z1+2].getParticleAt(0)->getF().operator[](0)==0);
    EXPECT_TRUE(c.Cells[2*y1*z1+z1+2].getParticleAt(0)->getF().operator[](1)==120);
    EXPECT_TRUE(c.Cells[2*y1*z1+z1+2].getParticleAt(0)->getF().operator[](2)==0);
}
TEST(hallocellsdecisions,y2reflective) {
    double sigma = 1;
    double epsilon = 5;
    double m = 1;
    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    double position[] = {0, 1, 0};
    double v[] = {0, 0, 0};
    Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), m, 0, epsilon, sigma);
    empty.push_back(a);
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    c.reassignMesh();
    ASSERT_TRUE(c.particles.size() == 1);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.y2boundary=2;
    c.calculateF();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    ASSERT_TRUE(c.particles.size() == 1);
    EXPECT_TRUE(c.Cells[2*y1*z1+3*z1+2].getElementNumber()==1);
    EXPECT_TRUE(c.Cells[2*y1*z1+3*z1+2].getParticleAt(0)->getF().operator[](0)==0);
    EXPECT_TRUE(c.Cells[2*y1*z1+3*z1+2].getParticleAt(0)->getF().operator[](1)==-120);
    EXPECT_TRUE(c.Cells[2*y1*z1+3*z1+2].getParticleAt(0)->getF().operator[](2)==0);
}
TEST(hallocellsdecisions,z1reflective) {
    double sigma = 1;
    double epsilon = 5;
    double m = 1;
    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    double position[] = {0, 0, -1};
    double v[] = {0, 0, 0};
    Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), m, 0, epsilon, sigma);
    empty.push_back(a);
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    c.reassignMesh();
    ASSERT_TRUE(c.particles.size() == 1);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.z1boundary=2;
    c.calculateF();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    ASSERT_TRUE(c.particles.size() == 1);
    EXPECT_TRUE(c.Cells[2*y1*z1+2*z1+1].getElementNumber()==1);
    EXPECT_TRUE(c.Cells[2*y1*z1+2*z1+1].getParticleAt(0)->getF().operator[](0)==0);
    EXPECT_TRUE(c.Cells[2*y1*z1+2*z1+1].getParticleAt(0)->getF().operator[](1)==0);
    EXPECT_TRUE(c.Cells[2*y1*z1+2*z1+1].getParticleAt(0)->getF().operator[](2)==120);
}
TEST(hallocellsdecisions,z2reflective) {
    double sigma = 1;
    double epsilon = 5;
    double m = 1;
    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    double position[] = {0, 0, 1};
    double v[] = {0, 0, 0};
    Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), m, 0, epsilon, sigma);
    empty.push_back(a);
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    c.reassignMesh();
    ASSERT_TRUE(c.particles.size() == 1);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.z2boundary=2;
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    c.calculateF();
    ASSERT_TRUE(c.particles.size() == 1);
    EXPECT_TRUE(c.Cells[2*y1*z1+2*z1+3].getElementNumber()==1);
    EXPECT_TRUE(c.Cells[2*y1*z1+2*z1+3].getParticleAt(0)->getF().operator[](0)==0);
    EXPECT_TRUE(c.Cells[2*y1*z1+2*z1+3].getParticleAt(0)->getF().operator[](1)==0);
    EXPECT_TRUE(c.Cells[2*y1*z1+2*z1+3].getParticleAt(0)->getF().operator[](2)==-120);
}
TEST(hallocellsdecisions,cornersofthecubeallreflective) {
    double sigma = 1;
    double epsilon = 5;
    double m = 1;
    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    double v[] = {0, 0, 0};
    for (int i = 1; i <= 2; i++) {
        for (int j = 1; j <= 2; j++) {
            for (int l = 1; l <= 2; l++) {
                double position[] = {pow(-1,i), pow(-1,j), pow(-1,l)};
                Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), m, 0, epsilon,
                                      sigma);
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    c.reassignMesh();
    ASSERT_TRUE(c.particles.size() == 8);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.x1boundary=2;
    c.x2boundary=2;
    c.y1boundary=2;
    c.y2boundary=2;
    c.z1boundary=2;
    c.z2boundary=2;
    c.calculateF();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    ASSERT_TRUE(c.particles.size() == 8);
    EXPECT_TRUE(c.Cells[y1*z1+3*z1+1].getElementNumber()==1);
    EXPECT_TRUE(c.Cells[y1*z1+3*z1+1].getParticleAt(0)->getF().operator[](0)==120);
    EXPECT_TRUE(c.Cells[y1*z1+3*z1+1].getParticleAt(0)->getF().operator[](1)==-120);
    EXPECT_TRUE(c.Cells[y1*z1+3*z1+1].getParticleAt(0)->getF().operator[](2)==120);
    EXPECT_TRUE(c.Cells[3*y1*z1+z1+1].getElementNumber()==1);
    EXPECT_TRUE(c.Cells[3*y1*z1+z1+1].getParticleAt(0)->getF().operator[](0)==-120);
    EXPECT_TRUE(c.Cells[3*y1*z1+z1+1].getParticleAt(0)->getF().operator[](1)==120);
    EXPECT_TRUE(c.Cells[3*y1*z1+z1+1].getParticleAt(0)->getF().operator[](2)==120);
    EXPECT_TRUE(c.Cells[y1*z1+z1+3].getElementNumber()==1);
    EXPECT_TRUE(c.Cells[y1*z1+z1+3].getParticleAt(0)->getF().operator[](0)==120);
    EXPECT_TRUE(c.Cells[y1*z1+z1+3].getParticleAt(0)->getF().operator[](1)==120);
    EXPECT_TRUE(c.Cells[y1*z1+z1+3].getParticleAt(0)->getF().operator[](2)==-120);
    EXPECT_TRUE(c.Cells[y1*z1+z1+1].getElementNumber()==1);
    EXPECT_TRUE(c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](0)==120);
    EXPECT_TRUE(c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](1)==120);
    EXPECT_TRUE(c.Cells[y1*z1+z1+1].getParticleAt(0)->getF().operator[](2)==120);
    EXPECT_TRUE(c.Cells[3*y1*z1+3*z1+1].getElementNumber()==1);
    EXPECT_TRUE(c.Cells[3*y1*z1+3*z1+1].getParticleAt(0)->getF().operator[](0)==-120);
    EXPECT_TRUE(c.Cells[3*y1*z1+3*z1+1].getParticleAt(0)->getF().operator[](1)==-120);
    EXPECT_TRUE(c.Cells[3*y1*z1+3*z1+1].getParticleAt(0)->getF().operator[](2)==120);
    EXPECT_TRUE(c.Cells[3*y1*z1+z1+3].getElementNumber()==1);
    EXPECT_TRUE(c.Cells[3*y1*z1+z1+3].getParticleAt(0)->getF().operator[](0)==-120);
    EXPECT_TRUE(c.Cells[3*y1*z1+z1+3].getParticleAt(0)->getF().operator[](1)==120);
    EXPECT_TRUE(c.Cells[3*y1*z1+z1+3].getParticleAt(0)->getF().operator[](2)==-120);
    EXPECT_TRUE(c.Cells[y1*z1+3*z1+3].getElementNumber()==1);
    EXPECT_TRUE(c.Cells[y1*z1+3*z1+3].getParticleAt(0)->getF().operator[](0)==120);
    EXPECT_TRUE(c.Cells[y1*z1+3*z1+3].getParticleAt(0)->getF().operator[](1)==-120);
    EXPECT_TRUE(c.Cells[y1*z1+3*z1+3].getParticleAt(0)->getF().operator[](2)==-120);
    EXPECT_TRUE(c.Cells[3*y1*z1+3*z1+3].getElementNumber()==1);
    EXPECT_TRUE(c.Cells[3*y1*z1+3*z1+3].getParticleAt(0)->getF().operator[](0)==-120);
    EXPECT_TRUE(c.Cells[3*y1*z1+3*z1+3].getParticleAt(0)->getF().operator[](1)==-120);
    EXPECT_TRUE(c.Cells[3*y1*z1+3*z1+3].getParticleAt(0)->getF().operator[](2)==-120);
}
TEST(hallocellsdecisions,x1periodic) {

    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    for (int i = -2; i <= 0; i++) {
        for (int j = -1; j <= 1; j++) {
            for (int l = -1; l <= 1; l++) {
                double c = i;
                double d = j;
                double e = l;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), 1, 0, 0, 0);
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size() == 27);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    for (int i = 0; i <= 2; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    c.x1boundary = 3;
    c.halloCellsDecisions();
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    for (int j=1;j<=3;j++){
        for(int l=1;l<=3;l++)
        {
            EXPECT_TRUE(c.Cells[j*z1+l].getElementNumber()==0);
        }
    }
    for (int j=1;j<=3;j++){
        for(int l=1;l<=3;l++)
        {
            EXPECT_TRUE(c.Cells[3*y1*z1+j*z1+l].getParticleAt(0)->getX().operator[](0)==1);
        }
    }
}
TEST(hallocellsdecisions,x2periodic) {

    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    for (int i = 0; i <= 2; i++) {
        for (int j = -1; j <= 1; j++) {
            for (int l = -1; l <= 1; l++) {
                double c = i;
                double d = j;
                double e = l;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), 1, 0, 0, 0);
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size() == 27);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    for (int i = 2; i <= 4; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    c.x2boundary = 3;
    c.halloCellsDecisions();
    ASSERT_TRUE(c.particles.size() == 27);
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    for (int j=1;j<=3;j++){
        for(int l=1;l<=3;l++)
        {
            EXPECT_TRUE(c.Cells[4*y1*z1+j*z1+l].getElementNumber()==0);
        }
    }
    for (int j=1;j<=3;j++){
        for(int l=1;l<=3;l++)
        {
            EXPECT_TRUE(c.Cells[1*y1*z1+j*z1+l].getParticleAt(0)->getX().operator[](0)==-1);
        }
    }
}
TEST(hallocellsdecisions,y1periodic) {

    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    for (int i = -1; i <= 1; i++) {
        for (int j = -2; j <= 0; j++) {
            for (int l = -1; l <= 1; l++) {
                double c = i;
                double d = j;
                double e = l;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), 1, 0, 0, 0);
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size() == 27);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    for (int i = 1; i <= 3; i++) {
        for (int j = 0; j <= 2; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    c.y1boundary = 3;
    c.halloCellsDecisions();
    ASSERT_TRUE(c.particles.size() == 27);
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    for (int i=1;i<=3;i++){
        for(int l=1;l<=3;l++)
        {
            EXPECT_TRUE(c.Cells[i*y1*z1+0*z1+l].getElementNumber()==0);
        }
    }
    for (int i=1;i<=3;i++){
        for(int l=1;l<=3;l++)
        {
            EXPECT_TRUE(c.Cells[i*y1*z1+3*z1+l].getParticleAt(0)->getX().operator[](1)==1);
        }
    }
}
TEST(hallocellsdecisions,y2periodic) {

    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    for (int i = -1; i <= 1; i++) {
        for (int j = 0; j <= 2; j++) {
            for (int l = -1; l <= 1; l++) {
                double c = i;
                double d = j;
                double e = l;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), 1, 0, 0, 0);
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size() == 27);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    for (int i = 1; i <= 3; i++) {
        for (int j = 2; j <= 4; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    c.y2boundary = 3;
    c.halloCellsDecisions();
    ASSERT_TRUE(c.particles.size() == 27);
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    for (int i=1;i<=3;i++){
        for(int l=1;l<=3;l++)
        {
            EXPECT_TRUE(c.Cells[i*y1*z1+4*z1+l].getElementNumber()==0);
        }
    }
    for (int i=1;i<=3;i++){
        for(int l=1;l<=3;l++)
        {
            EXPECT_TRUE(c.Cells[i*y1*z1+z1+l].getParticleAt(0)->getX().operator[](1)==-1);
        }
    }
}
TEST(hallocellsdecisions,z1periodic) {

    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            for (int l = -2; l <= 0; l++) {
                double c = i;
                double d = j;
                double e = l;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), 1, 0, 0, 0);
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size() == 27);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 0; l <= 2; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    c.z1boundary = 3;
    c.halloCellsDecisions();
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    for (int i=1;i<=3;i++){
        for(int j=1;j<=3;j++)
        {
            EXPECT_TRUE(c.Cells[i*y1*z1+j*z1].getElementNumber()==0);
        }
    }
    for (int i=1;i<=3;i++){
        for(int j=1;j<=3;j++)
        {
            EXPECT_TRUE(c.Cells[i*y1*z1+j*z1+3].getParticleAt(0)->getX().operator[](2)==1);
        }
    }
}
TEST(hallocellsdecisions,z2periodic) {

    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            for (int l = 0; l <= 2; l++) {
                double c = i;
                double d = j;
                double e = l;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), 1, 0, 0, 0);
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size() == 27);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 2; l <= 4; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    c.z2boundary = 3;
    c.halloCellsDecisions();
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    for (int i=1;i<=3;i++){
        for(int j=1;j<=3;j++)
        {
            EXPECT_TRUE(c.Cells[i*y1*z1+j*z1+4].getElementNumber()==0);
        }
    }
    for (int i=1;i<=3;i++){
        for(int j=1;j<=3;j++)
        {
            EXPECT_TRUE(c.Cells[i*y1*z1+j*z1+1].getParticleAt(0)->getX().operator[](2)==-1);
        }
    }
}
TEST(hallocellsdecisions,allperiodic) {

    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    for (int i = -2; i <= 2; i=i+2) {
        for (int j = -2; j <= 2; j=j+2) {
            for (int l = -2; l <= 2; l=l+2) {
                    double c = i;
                    double d = j;
                    double e = l;
                    double position[] = {c, d, e};
                    double v[] = {0, 0, 0};
                    Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), 1, 0, 0, 0);
                    empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size() == 27);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    for (int i = 0; i <= 4; i=i+2) {
        for (int j = 0; j <= 4; j=j+2) {
            for (int l = 0; l <= 4; l=l+2) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    c.x1boundary = 3;
    c.x2boundary = 3;
    c.y1boundary = 3;
    c.y2boundary = 3;
    c.z1boundary = 3;
    c.z2boundary = 3;
    c.halloCellsDecisions();
    ASSERT_TRUE(c.particles.size() == 27);
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    for (int i=1;i<=3;i++){
        for(int j=1;j<=3;j++)
        {
            EXPECT_TRUE(c.Cells[i*y1*z1+j*z1+4].getElementNumber()==0);
            EXPECT_TRUE(c.Cells[i*y1*z1+j*z1+0].getElementNumber()==0);
            EXPECT_TRUE(c.Cells[i*y1*z1+0*z1+j].getElementNumber()==0);
            EXPECT_TRUE(c.Cells[i*y1*z1+4*z1+j].getElementNumber()==0);
            EXPECT_TRUE(c.Cells[0*y1*z1+i*z1+j].getElementNumber()==0);
            EXPECT_TRUE(c.Cells[4*y1*z1+i*z1+j].getElementNumber()==0);
        }
    }
    for (int i=1;i<=3;i++){
        for(int j=1;j<=3;j++)
        {
            EXPECT_TRUE(c.Cells[i*y1*z1+j*z1+1].getParticleAt(0)->getX().operator[](2)==-1);
            EXPECT_TRUE(c.Cells[i*y1*z1+j*z1+3].getParticleAt(0)->getX().operator[](2)==1);
            EXPECT_TRUE(c.Cells[i*y1*z1+1*z1+j].getParticleAt(0)->getX().operator[](1)==-1);
            EXPECT_TRUE(c.Cells[i*y1*z1+3*z1+j].getParticleAt(0)->getX().operator[](1)==1);
            EXPECT_TRUE(c.Cells[1*y1*z1+i*z1+j].getParticleAt(0)->getX().operator[](0)==-1);
            EXPECT_TRUE(c.Cells[3*y1*z1+i*z1+j].getParticleAt(0)->getX().operator[](0)==1);
        }
    }
}
TEST(hallocellsdecisions,x1x2periodicforce) {

    double sigma = 1;
    double epsilon = 5;
    double m = 1;
    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    for (int i = -1; i <= 1; i=i+2) {
        for (int j = -1; j <= 1; j=j+1) {
            for (int l = -1; l <= 1; l=l+1) {
                double c = i;
                double d = j;
                double e = l;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), m, 0, epsilon, sigma);
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size() == 18);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    for(int i=1;i<=3;i=i+2) {
        for(int j=1;j<=3;j=j+1) {
            for(int l=1;l<=3;l=l+1) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    c.x1boundary = 3;
    c.x2boundary = 3;
    c.calculateF();
    for(int i=1;i<=3;i=i+2) {
        for(int j=1;j<=3;j=j+1) {
            for(int l=1;l<=3;l=l+1) {
                for(int k=0;k<c.Cells[i*y1*z1+j*z1+l].getElementNumber();k++){
                    EXPECT_EQ(c.Cells[i*y1*z1+j*z1+l].getParticleAt(k)->getF().operator[](0),(2-i)*(120));
                    if(c.Cells[i*y1*z1+j*z1+l].getParticleAt(k)->getF().operator[](0)!=(2-i)*(120))
                        std::cout<<i<<" "<<j<<" "<<l<<std::endl;
                }
            }
        }
    }

}
TEST(hallocellsdecisions,y1y2periodicforce) {

    double sigma = 1;
    double epsilon = 5;
    double m = 1;
    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    for (int i = -1; i <= 1; i=i+1) {
        for (int j = -1; j <= 1; j=j+2) {
            for (int l = -1; l <= 1; l=l+1) {
                double c = i;
                double d = j;
                double e = l;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), m, 0, epsilon, sigma);
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size() == 18);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    for(int i=1;i<=3;i=i+1) {
        for(int j=1;j<=3;j=j+2) {
            for(int l=1;l<=3;l=l+1) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    c.y1boundary = 3;
    c.y2boundary = 3;
    c.calculateF();
    for(int i=1;i<=3;i=i+1) {
        for(int j=1;j<=3;j=j+2) {
            for(int l=1;l<=3;l=l+1) {
                for(int k=0;k<c.Cells[i*y1*z1+j*z1+l].getElementNumber();k++){
                    EXPECT_EQ(c.Cells[i*y1*z1+j*z1+l].getParticleAt(k)->getF().operator[](1),(2-j)*(120));
                    if(c.Cells[i*y1*z1+j*z1+l].getParticleAt(k)->getF().operator[](1)!=(2-j)*(120))
                        std::cout<<i<<" "<<j<<" "<<2<<std::endl;
                }
            }
        }
    }

}
TEST(hallocellsdecisions,z1z2periodicforce) {

    double sigma = 1;
    double epsilon = 5;
    double m = 1;
    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    for (int i = -1; i <= 1; i=i+1) {
        for (int j = -1; j <= 1; j=j+1) {
            for (int l = -1; l <= 1; l=l+2) {
                double c = i;
                double d = j;
                double e = l;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), m, 0, epsilon, sigma);
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size() == 18);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    for(int i=1;i<=3;i=i+1) {
        for(int j=1;j<=3;j=j+1) {
            for(int l=1;l<=3;l=l+2) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    c.z1boundary = 3;
    c.z2boundary = 3;
    c.calculateF();
    for(int i=1;i<=3;i=i+1) {
        for(int j=1;j<=3;j=j+1) {
            for(int l=1;l<=3;l=l+2) {
                for(int k=0;k<c.Cells[i*y1*z1+j*z1+l].getElementNumber();k++){
                    EXPECT_EQ(c.Cells[i*y1*z1+j*z1+l].getParticleAt(k)->getF().operator[](2),(2-l)*(120));
                }
            }
        }
    }

}
TEST(hallocellsdecisions,allperiodicforce) {

    double sigma = 1;
    double epsilon = 5;
    double m = 1;
    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    for (int i = -1; i <= 1; i=i+2) {
        for (int j = -1; j <= 1; j=j+2) {
            for (int l = -1; l <= 1; l=l+2) {
                double c = i;
                double d = j;
                double e = l;
                double position[] = {c, d, e};
                double v[] = {0, 0, 0};
                Particle a = Particle(utils::Vector<double, 3>(position), utils::Vector<double, 3>(v), m, 0, epsilon, sigma);
                empty.push_back(a);
            }
        }
    }
    utils::Vector<double, 3> d = utils::Vector<double, 3>(x);
    Mesh c = Mesh(false,0,1, d, empty, 1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size() == 8);
    ASSERT_TRUE(c.xsize == 5);
    ASSERT_TRUE(c.ysize == 5);
    ASSERT_TRUE(c.zsize == 5);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    for(int i=1;i<=3;i=i+2) {
        for(int j=1;j<=3;j=j+2) {
            for(int l=1;l<=3;l=l+2) {
                ASSERT_TRUE(c.Cells[i*y1*z1+j*z1+l].getElementNumber() == 1);
            }
        }
    }
    c.x1boundary = 3;
    c.x2boundary = 3;
    c.y1boundary = 3;
    c.y2boundary = 3;
    c.z1boundary = 3;
    c.z2boundary = 3;
    c.calculateF();
    for(int i=1;i<=3;i=i+2) {
        for(int j=1;j<=3;j=j+2) {
            for(int l=1;l<=3;l=l+2) {
                for(int k=0;k<c.Cells[i*y1*z1+j*z1+l].getElementNumber();k++){
                    EXPECT_EQ(c.Cells[i*y1*z1+j*z1+l].getParticleAt(k)->getF().operator[](0),(2-i)*(120));
                    EXPECT_EQ(c.Cells[i*y1*z1+j*z1+l].getParticleAt(k)->getF().operator[](1),(2-j)*(120));
                    EXPECT_EQ(c.Cells[i*y1*z1+j*z1+l].getParticleAt(k)->getF().operator[](2),(2-l)*(120));
                    if(c.Cells[i*y1*z1+j*z1+l].getParticleAt(k)->getF().operator[](1)!=(2-j)*(120)||c.Cells[i*y1*z1+j*z1+l].getParticleAt(k)->getF().operator[](2)!=(2-l)*(120))
                        std::cout<<i<<" "<<j<<" "<<l<<std::endl;
                }
            }
        }
    }

}
TEST(testF, TestFnewton3) {
    double domain[] = {303.0, 180.0, 0};
    utils::Vector<double, 3> domainsize = utils::Vector<double, 3>(domain);
    double cutoff = 3.0;
    double d_t = 0.00005;
    int dimension = 2;
    double sigma = 1.2;
    double epsilon = 1.0;
    std::list<Particle> particles;
    double radius = 20.0;
    double h = 1.2;
    double factor = 0.1;
    double centre[] = {150.0 - domain[0] / 2, 150.0 - domain[1] / 2, 1.0 - domain[2] / 2};
    double v[] = {0.0, 0.0, 0.0};
    double mass = 1.0;
    utils::Vector<double, 3> velocity = utils::Vector<double, 3>(v);
    utils::Vector<double, 3> coordinates = utils::Vector<double, 3>(centre);
    double mesh_width = h;
    utils::Vector<double, 3> dimensions = dimensions.operator=(ceil((radius * 2) / h));
    for (int i = 0; i < dimensions[0]; i++) {
        double a1 = coordinates[0] - radius + mesh_width * i;
        for (int j = 0; j < dimensions[1]; j++) {
            double a2 = coordinates[1] - radius + mesh_width * j;
            if (dimension == 3)
                for (int k = 0; k < dimensions[2]; k++) {
                    double a[3] = {a1, a2, coordinates[2] - radius + mesh_width * k};
                    utils::Vector<double, 3> c2(a);
                    utils::Vector<double, 3> c3(a);
                    c3.operator=(c3.operator-(coordinates));
                    if (c3.L2Norm() <= radius) {
                        Particle p(c2, velocity, mass, 0, epsilon, sigma);
                        MaxwellBoltzmannDistribution(p, factor, 2);
                        particles.push_back(p);
                    }
                }
            else {
                double a[3] = {a1, a2, 0};
                utils::Vector<double, 3> c2(a);
                utils::Vector<double, 3> c3(a);
                c3.operator=(c3.operator-(coordinates));
                ///@brief point is actually generated only if the point is in the radius of the intended circle/sphere
                if (c3.L2Norm() <= radius) {
                    Particle p(c2, velocity, (i%2)*mass, 0, epsilon, sigma);
                    MaxwellBoltzmannDistribution(p, factor * sqrt(8 / 3.14), 2);
                    particles.push_back(p);
                }
            }
        }
    }
    std::cout << particles.size() << std::endl;
    Mesh CellMesh = Mesh(false,0,cutoff, domainsize, particles, d_t, -12.44, 2, 2, 2, 2, 0, 0);
    ASSERT_TRUE(CellMesh.particles.size() == particles.size());
    ASSERT_TRUE(CellMesh.xsize == 103);
    ASSERT_TRUE(CellMesh.ysize == 62);
    ASSERT_TRUE(CellMesh.zsize == 3);
    CellMesh.reassignMesh();
    int x1=CellMesh.xsize;
    int y1=CellMesh.ysize;
    int z1=CellMesh.zsize;
    for (int i = 0; i <50;i++){
        std::cout<<i<<std::endl;
        CellMesh.calculateX();
        CellMesh.calculateF();
        CellMesh.calculateV();
    utils::Vector<double, 3> sum = utils::Vector<double, 3>(0.0);
    int vis = 0;
    for (int i = 0; i < CellMesh.xsize; i++) {
        for (int j = 0; j < CellMesh.ysize; j++) {
            for (int l = 0; l < CellMesh.zsize; l++) {

                for (int k = 0; k < CellMesh.Cells[i*y1*z1+j*z1+l].getElementNumber(); k++) {
                    sum.operator=(sum.operator+(CellMesh.Cells[i*y1*z1+j*z1+l].getParticleAt(k)->getF()));
                    vis = vis + 1;
                }
            }
        }
    }
    EXPECT_TRUE(sum.operator[](2) == 0);
    EXPECT_NEAR(sum.operator[](0), 0, 1 * pow(10, -10));
    EXPECT_NEAR(vis *(3/2)* -12.44, sum.operator[](1), 1 * pow(10, -10));
}

}
TEST(testv,testforce2particlesinsamecellatdistancegreaterthan6thsquarerootofsigmaandcalculatev)
{

    double sigma=1;
    double epsilon=5;
    double m=1;
    std::list<Particle> empty;
    double x[]={2,2,0};
    utils::Vector<double,3> d=utils::Vector<double ,3>(x);
    double positionp[] = {0.5, 0.5, 0};
    double positionn[] = {-0.5, -0.5, 0};
    double v[] = {0, 0, 0};
    Particle a = Particle(utils::Vector<double ,3>(positionp),utils::Vector<double ,3>(v), m, 0, epsilon, sigma );
    empty.push_front(a);
    Particle b = Particle(utils::Vector<double ,3>(positionn),utils::Vector<double ,3>(v), m, 0, epsilon, sigma );
    empty.push_front(b);
    Mesh c=Mesh(false,0,2.0,d,empty,1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size()==2);
    ASSERT_TRUE(c.xsize==3);
    ASSERT_TRUE(c.ysize==3);
    ASSERT_TRUE(c.zsize==3);
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    c.reassignMesh();
    ASSERT_TRUE(c.Cells[1*y1*z1+1*z1+1].getElementNumber()==2);
    ASSERT_TRUE((c.Cells[1*y1*z1+1*z1+1].getParticleAt(0)->getX().operator-(c.Cells[1*y1*z1+1*z1+1].getParticleAt(1)->getX())).L2Norm()>(pow(2,1/6)*sigma));
    c.calculateF();
    ASSERT_TRUE(c.Cells[1*y1*z1+1*z1+1].getParticleAt(0)->getF().operator==(c.Cells[1*y1*z1+1*z1+1].getParticleAt(1)->getF().operator*(-1)));
    EXPECT_TRUE(c.Cells[1*y1*z1+1*z1+1].getParticleAt(0)->getF().operator[](0)>0&&c.Cells[1*y1*z1+1*z1+1].getParticleAt(0)->getF().operator[](1)>0&&c.Cells[1*y1*z1+1*z1+1].getParticleAt(0)->getF().operator[](2)==0);
    EXPECT_TRUE(c.Cells[1*y1*z1+1*z1+1].getParticleAt(1)->getF().operator[](0)<0&&c.Cells[1*y1*z1+1*z1+1].getParticleAt(1)->getF().operator[](1)<0&&c.Cells[1*y1*z1+1*z1+1].getParticleAt(1)->getF().operator[](2)==0);
    c.calculateV();
    ASSERT_TRUE(c.Cells[1*y1*z1+1*z1+1].getParticleAt(0)->getV().operator==(c.Cells[1*y1*z1+1*z1+1].getParticleAt(1)->getV().operator*(-1)));
    EXPECT_TRUE(c.Cells[1*y1*z1+1*z1+1].getParticleAt(0)->getV().operator[](0)>0&&c.Cells[1*y1*z1+1*z1+1].getParticleAt(0)->getV().operator[](1)>0&&c.Cells[1*y1*z1+1*z1+1].getParticleAt(0)->getV().operator[](2)==0);
    EXPECT_TRUE(c.Cells[1*y1*z1+1*z1+1].getParticleAt(1)->getV().operator[](0)<0&&c.Cells[1*y1*z1+1*z1+1].getParticleAt(1)->getV().operator[](1)<0&&c.Cells[1*y1*z1+1*z1+1].getParticleAt(1)->getV().operator[](2)==0);
}
TEST(CalculateX,Celltransition) {

    double sigma = 1;
    double epsilon = 5;
    double m = 1;
    std::list<Particle> empty;
    double x[] = {3, 3, 3};
    utils::Vector<double,3> d=utils::Vector<double ,3>(x);
    double positionp[] = {1, 1, 1};
    double v[] = {-1, -1, -1};
    Particle a = Particle(utils::Vector<double ,3>(positionp),utils::Vector<double ,3>(v), m, 0, epsilon, sigma );
    empty.push_front(a);
    Mesh c=Mesh(false,0,1.0,d,empty,1,0,0,0,0,0,0,0);
    ASSERT_TRUE(c.particles.size()==1);
    ASSERT_TRUE(c.xsize==5);
    ASSERT_TRUE(c.ysize==5);
    ASSERT_TRUE(c.zsize==5);
    c.reassignMesh();
    int x1=c.xsize;
    int y1=c.ysize;
    int z1=c.zsize;
    ASSERT_TRUE(c.Cells[3*y1*z1+3*z1+3].getElementNumber()==1);
    c.calculateX();
    EXPECT_TRUE(c.Cells[2*y1*z1+2*z1+2].getElementNumber()==1);
    EXPECT_TRUE(c.Cells[2*y1*z1+2*z1+2].getParticleAt(0)->getX().operator==(utils::Vector<double,3>(0.0)));

}
int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
