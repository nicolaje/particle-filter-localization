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
