cmake_minimum_required(VERSION 2.8)

#SET(CMAKE_BUILD_TYPE distribution)
SET(CMAKE_BUILD_TYPE debug)
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-g -O3")
SET(CMAKE_C_FLAGS_DISTRIBUTION "-g -O3")
SET(CMAKE_CXX_FLAGS_DEBUG "-g")
SET(CMAKE_C_FLAGS_DEBUG "-g")

set(CMAKE_VERBOSE_MAKEFILE ON)

project(bp-vision)
#LINK   = -lm 

include_directories(.)
add_executable(restore restore.cpp)
add_executable(stereo stereo.cpp)
add_executable(depth depth.cpp)
add_executable(noise noise.cpp)

#add_library(test STATIC test.c)
#target_link_libraries(myapp test)
#install(TARGETS test DESTINATION lib)

#install(TARGETS myapp DESTINATION bin)
