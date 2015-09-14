CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

ADD_DEFINITIONS(
  -std=c++11
)

FIND_PACKAGE(octomap REQUIRED)

INCLUDE_DIRECTORIES("Eigen/")

project(lib_particle_filter)

LIST(APPEND CMAKE_INCLUDE_PATH ${CMAKE_CURRENT_LIST_DIR}/Eigen)

ADD_LIBRARY(particle_filter SHARED src/ParticleFilter.cpp)
TARGET_LINK_LIBRARIES(particle_filter octomap)


ADD_EXECUTABLE(example src/example.cpp src/simulator.cpp src/ParticleFilter.cpp src/vibes.cpp src/simulator.h src/ParticleFilter.h src/vibes.h)
TARGET_LINK_LIBRARIES(example particle_filter octomap)