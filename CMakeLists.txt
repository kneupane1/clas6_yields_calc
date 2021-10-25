cmake_minimum_required(VERSION 3.12)
project(PhysicsEventClasses)

set(CMAKE_CXX_STANDARD 17)

include_directories(src)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO) #RooFit RooFitCore RooStats RootAuth)
include(${ROOT_USE_FILE})

add_executable(PhysicsEventClasses
        src/branches.cpp
        src/branches.hpp
        src/clipp.h
        src/color.hpp
        src/constants.hpp
        src/histogram.cpp
        src/histogram.hpp
        src/main.cpp
        src/main.hpp
        src/physics.cpp
        src/physics.hpp
        src/reaction.cpp
        src/reaction.hpp
        # clas6_yields.cpp
        # clas12_yields.cpp
        )

target_link_libraries(PhysicsEventClasses ${ROOT_LIBRARIES})