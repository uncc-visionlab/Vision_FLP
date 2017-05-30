cmake_minimum_required(VERSION 2.8)

SET(CMAKE_BUILD_TYPE distribution)
#SET(CMAKE_BUILD_TYPE debug)
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-g -W -Wall -O3 -DUSE_64_BIT_PTR_CAST")
SET(CMAKE_C_FLAGS_DISTRIBUTION "-g -W -Wall -O3")
SET(CMAKE_CXX_FLAGS_DEBUG "-g -W -Wall -DUSE_64_BIT_PTR_CAST")
SET(CMAKE_C_FLAGS_DEBUG "-g -W -Wall")

project(MRF2.2)

include_directories(.)

#SRC =  mrf.cpp ICM.cpp GCoptimization.cpp graph.cpp maxflow.cpp \
#       MaxProdBP.cpp LinkedBlockList.cpp regions-maxprod.cpp \
#       TRW-S.cpp BP-S.cpp

add_library(MRF STATIC mrf.cpp ICM.cpp GCoptimization.cpp graph.cpp maxflow.cpp
       MaxProdBP.cpp LinkedBlockList.cpp regions-maxprod.cpp
       TRW-S.cpp BP-S.cpp)

add_executable(example-mrf example.cpp)
target_link_libraries(example-mrf MRF)

#mrf.o: mrf.h
#ICM.o: ICM.h mrf.h LinkedBlockList.h
#GCoptimization.o: energy.h graph.h block.h mrf.h GCoptimization.h
#GCoptimization.o: LinkedBlockList.h
#graph.o: graph.h block.h mrf.h
#maxflow.o: graph.h block.h mrf.h
#MaxProdBP.o: MaxProdBP.h mrf.h LinkedBlockList.h regions-new.h
#LinkedBlockList.o: LinkedBlockList.h
#regions-maxprod.o: MaxProdBP.h mrf.h LinkedBlockList.h regions-new.h
#TRW-S.o: TRW-S.h mrf.h typeTruncatedQuadratic2D.h
#BP-S.o: BP-S.h mrf.h typeTruncatedQuadratic2D.h
