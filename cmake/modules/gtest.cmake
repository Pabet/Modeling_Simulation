enable_testing()
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/GoogleTestMaster/googletest
        ${CMAKE_CURRENT_SOURCE_DIR}/GoogleTestMaster
        EXCLUDE_FROM_ALL)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
add_executable(Calculatortest src/Calculator.cpp src/Calculator.h src/tests/CalculatorTest.cpp src/tests/CalculatorTest.h src/Particle.cpp src/Particle.h  )
target_include_directories(Calculatortest
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/libxsd
        )
target_link_libraries(Calculatortest PUBLIC gtest gtest_main)
add_test(Calculatortest Calculatortest)
add_executable(ParticleContainerTest src/particleContainer.cpp src/particleContainer.h src/tests/particleContainerTest.cpp  src/tests/particleContainerTest.h src/Particle.cpp src/Particle.h  )
target_include_directories(ParticleContainerTest
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/libxsd
        )
target_link_libraries(ParticleContainerTest PUBLIC gtest gtest_main)
add_test(ParticleContainerTest ParticleContainerTest)
add_executable(MeshTest src/Mesh.cpp src/Mesh.h src/tests/MeshTest.cpp  src/tests/MeshTest.h src/Particle.cpp src/Particle.h src/particleContainer.cpp src/particleContainer.h  src/Calculator.cpp src/Calculator.h src/MaxwellBoltzmannDistribution.h src/MaxwellBoltzmannDistribution.cpp)
target_include_directories(MeshTest
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/libxsd
        )
target_link_libraries(MeshTest PUBLIC gtest gtest_main)
add_test(MeshTest MeshTest)