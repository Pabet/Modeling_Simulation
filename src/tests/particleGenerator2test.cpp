//
// Created by Victor Stroescu on 05/12/2018.
//

#include "particleGenerator2test.h"
#include "gtest/gtest.h"
#include "particleGenerator2.h"
#include "utils/Vector.h"

TEST(CalculateForce, nootherparticles)
{
    utils::Vector<double,3> d=utils::Vector<double ,3>({120.0});
    utils::Vector<double,3> f=utils::Vector<double ,3>({0.0});
    Particle p=Particle(0);
    Particle n[20];
    n[0]=p;
    particleGenerator2 a= particleGenerator2(true, d, 10);
    particleGenerator2::molecule s={f,n};


}
int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
