cmake_minimum_required(VERSION 3.5)

#logging switch (ON/OFF) if disabled no logging library will be
#included and all loggin statements will be disabled before compilation
#LOGGING OFF to use programm on LRZ
set(LOGGING On)

# define project name, version, and language
project(PSEMolDyn_GroupG VERSION 0.0.1 LANGUAGES CXX)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# set Release as the default build type if it is not yet set.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# collect all cpp files
file(GLOB_RECURSE MY_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
        # header don't need to be included but this might be necessary for some IDEs
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
        )


# create make target
add_executable(MolSim src/outputWriter/vtk-unstructured.cpp src/outputWriter/XYZWriter.cpp src/Calculator.cpp src/MolSim.cpp src/MaxwellBoltzmannDistribution.cpp src/particleGenerator.cpp src/Particle.cpp src/FileReader.cpp src/Input.cpp src/Input.h src/outputWriter/pvdWriter.cpp src/outputWriter/VTKWriter.cpp src/outputWriter/pvdWriter.h src/shapes.hxx src/shapes.cxx src/particleGenerator2.cpp src/SettingsContainer.cpp src/SettingsContainer.h src/thermostat.cpp src/thermostat.h src/outputWriter/checkpoint.cpp src/outputWriter/checkpoint.h src/Mesh.cpp src/Mesh.h src/particleContainer.cpp src/particleContainer.h)


# set cxx standard. You may raise this if you want.
target_compile_features(MolSim
        PRIVATE
        cxx_std_11
        )

target_include_directories(MolSim
        PUBLIC
        SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/libxsd
        PRIVATE
        SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/src
        )

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(LOGGING)
    target_compile_definitions(MolSim
            PUBLIC
            logOn
            )
    target_link_libraries(MolSim
        PUBLIC
        xerces-c
        log4cxx
    )
else()
    target_link_libraries(MolSim
        PUBLIC
        xerces-c
)
endif()

# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
target_compile_options(MolSim
        PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-Wall>
        # disable some spam warnings for icpc...
        $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
        )

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
include(gtest)
include(doxygen)
