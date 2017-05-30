cmake_minimum_required(VERSION 2.8)

SET(CMAKE_BUILD_TYPE distribution)
#SET(CMAKE_BUILD_TYPE debug)
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-g -W -Wall -O3")
SET(CMAKE_C_FLAGS_DISTRIBUTION "-g -W -Wall -O3")
#SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-g -W -Wall -O3 -isystem")
#SET(CMAKE_C_FLAGS_DISTRIBUTION "-g -W -Wall -O3 -isystem")
SET(CMAKE_CXX_FLAGS_DEBUG "-g")
SET(CMAKE_C_FLAGS_DEBUG "-g")

project(TRW_S-v1.3)

include_directories(.)

add_library(MRF_TRW STATIC minimize.cpp MRFEnergy.cpp ordering.cpp treeProbabilities.cpp)

#add_executable(example-trw example.cpp)
#target_link_libraries(example-trw MRF_TRW)

