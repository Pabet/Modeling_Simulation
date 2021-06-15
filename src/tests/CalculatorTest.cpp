//
// Created by Victor Stroescu on 15/11/2018.
//

#include "CalculatorTest.h"
#include "gtest/gtest.h"
#include "Calculator.h"
#include "utils/Vector.h"

TEST(velocity, nootherforces){
    Calculator s;
    utils::Vector<double, 3> a= utils::Vector<double,3>(1.0);
    Particle T1=Particle(utils::Vector<double,3>(0.0),utils::Vector<double,3>(1.0),5,0);
    std::list<Particle> particles;
    particles.push_back(T1);
    s.calculateV(particles,6);
    ASSERT_FALSE(particles.size()==0);

    EXPECT_TRUE(a.operator==(particles.front().getV()));

}
TEST(force,nootherparticles){
    Calculator s;
    Particle T1=Particle(utils::Vector<double,3>(0.0),utils::Vector<double,3>(1.0),5,0);
    std::list<Particle> particles;
    particles.push_back(T1);
    s.calculateF(particles,0,0,0);
    utils::Vector<double, 3> a= utils::Vector<double,3>(0.0);
    ASSERT_FALSE(particles.size()==0);
    EXPECT_TRUE(a.operator==(particles.front().getF()));
}
TEST(position,simplevelocity)
{
    Calculator s;
    utils::Vector<double, 3> a= utils::Vector<double,3>(1.0);
    Particle T1=Particle(utils::Vector<double,3>(0.0),utils::Vector<double,3>(1.0),5,0);
    std::list<Particle> particles;
    particles.push_back(T1);
    s.calculateX(particles,1.0);
    ASSERT_FALSE(particles.size()==0);
    EXPECT_TRUE(a.operator==(particles.front().getX()));

}
int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}