//
// Created by Victor Stroescu on 12/01/2019.
//

#include "particleContainerTest.h"
#include "particleContainer.h"
#include "gtest/gtest.h"
#include "Particle.h"
#include "utils/Vector.h"
TEST(addtest,addandget1particletoempty)
{
    particleContainer a=particleContainer();
    double x []={0,1,2};
    Particle* s=new Particle[1];
    s[0]=Particle(utils::Vector<double,3>(x),utils::Vector<double,3>(5),2,0,0,0);
    a.addParticle(s);
    ASSERT_FALSE(a.getElementNumber()!=1);
    EXPECT_TRUE(s->operator==(*a.getParticleAt(0)));
}
TEST(addtest,addandget1particletoemptycomparedtocopy)
{
    particleContainer a=particleContainer();
    double x []={0,1,2};

    Particle* s=new Particle[1];
    Particle* scopy=new Particle[1];
    s[0]=Particle(utils::Vector<double,3>(x),utils::Vector<double,3>(5),2,0,0,0);
    scopy[0]=Particle(utils::Vector<double,3>(x),utils::Vector<double,3>(5),2,0,0,0);
    a.addParticle(s);
    ASSERT_FALSE(a.getElementNumber()!=1);
    EXPECT_TRUE(scopy->operator==(*a.getParticleAt(0)));
}
TEST(addtest,addandget2particlefromempty)
{
    particleContainer a=particleContainer();
    ASSERT_TRUE(a.getSize()==2);
    double x1 []={0,1,2};
    Particle* s1=new Particle[1];
    s1[0]=Particle(utils::Vector<double,3>(x1),utils::Vector<double,3>(5),2,0,0,0);
    double x2[]={6,7,8};
    Particle* s2=new Particle[1];
    s2[0]=Particle(utils::Vector<double,3>(x2),utils::Vector<double,3>(7),2,0,0,0);
    a.addParticle(s1);
    a.addParticle(s2);
    ASSERT_FALSE(a.getElementNumber()!=2);
    EXPECT_TRUE(s1->operator==(*a.getParticleAt(0)));
    EXPECT_TRUE(s2->operator==(*a.getParticleAt(1)));
    EXPECT_TRUE(a.getSize()==4);
}
TEST(removetest,add16elementsthenremovethe6th)
{
    particleContainer a=particleContainer();
    for(int i=1;i<=16;i++)
    {
        double k=(double) i;
        double x1[]={k,2*k,3*k};
        Particle* s1=new Particle[1];
        s1[0]=Particle(utils::Vector<double,3>(x1),utils::Vector<double,3>(4*k),5*k,1,0,0);
        a.addParticle(s1);
    }
    ASSERT_TRUE(a.getElementNumber()==16);
    ASSERT_TRUE(a.getSize()==32);
    for(int i=1;i<=a.getElementNumber();i++)
    {
        double k=(double) i;
        double x1[]={k,2*k,3*k};
        Particle s1=Particle(utils::Vector<double,3>(x1),utils::Vector<double,3>(4*k),5*k,1,0,0);
        ASSERT_TRUE(a.getParticleAt(i-1)->operator==(s1));
    }
    a.removeparticleAt(6);
    for(int i=1;i<=a.getElementNumber();i++)
    {
        double k=(double) i;
        if(i>=7)
        {
            k++;
        }
        double x1[]={k,2*k,3*k};
        Particle s1=Particle(utils::Vector<double,3>(x1),utils::Vector<double,3>(4*k),5*k,1,0,0);
        ASSERT_TRUE(a.getParticleAt(i-1)->operator==(s1));
    }
    ASSERT_TRUE(a.getSize()==16);
    ASSERT_TRUE(a.getElementNumber()==15);
}
TEST(addtest,add16elements)
{
    particleContainer a=particleContainer();
    for(int i=1;i<=16;i++)
    {
        double k=(double) i;
        double x1[]={k,2*k,3*k};
        Particle* s1=new Particle[1];
        s1[0]=Particle(utils::Vector<double,3>(x1),utils::Vector<double,3>(4*k),5*k,1,0,0);
        a.addParticle(s1);
    }
    ASSERT_TRUE(a.getElementNumber()==16);
    EXPECT_TRUE(a.getSize()==32);
    for(int i=1;i<=a.getElementNumber();i++)
    {
        double k=(double) i;
        double x1[]={k,2*k,3*k};
        Particle s1=Particle(utils::Vector<double,3>(x1),utils::Vector<double,3>(4*k),5*k,1,0,0);
        EXPECT_TRUE(a.getParticleAt(i-1)->operator==(s1));
    }
}
TEST(gettest,getunaddedparticlefromempty)
{
    particleContainer a=particleContainer();
    ASSERT_TRUE(a.getSize()==2);
    ASSERT_TRUE(a.getElementNumber()==0);
    EXPECT_THROW(a.getParticleAt(0),std::runtime_error);
}
int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
