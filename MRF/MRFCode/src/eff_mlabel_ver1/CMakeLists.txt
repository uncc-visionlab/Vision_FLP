cmake_minimum_required(VERSION 2.8)

#SET(CMAKE_BUILD_TYPE distribution)
SET(CMAKE_BUILD_TYPE debug)
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-g -O3 -DNDEBUG")# -Wformat=")
SET(CMAKE_C_FLAGS_DISTRIBUTION "-g -O3 -DNDEBUG")# -Wformat=")
SET(CMAKE_CXX_FLAGS_DEBUG "-g")
SET(CMAKE_C_FLAGS_DEBUG "-g")

project(eff_mlabel_ver1)

include_directories(.)

add_executable(example-eff-mlabel main.cpp image.cpp)

#CCOMP = g++-4.0
#CFLAGS_DEBUG = -g -pg
#CFLAGS_RELEASE = -DNDEBUG -pg
#CFLAGS = $(CFLAGS_DEBUG)
#OBJ_DIR = ./.obj
#OUTPUT = ./minimize

#all: $(OUTPUT)

#$(OUTPUT): $(OBJ_DIR)/main.o $(OBJ_DIR)/image.o
#	$(CCOMP) $(CFLAGS) -lm -o $(OUTPUT) $(OBJ_DIR)/main.o $(OBJ_DIR)/image.o

#$(OBJ_DIR)/main.o: main.cpp image.h includes.h energy.h graph.h kovtun.h aexpand.h TRWBP.h typePotts.h MRFEnergy.h block.h
#	$(CCOMP) $(CFLAGS) -c -o $(OBJ_DIR)/main.o main.cpp

#$(OBJ_DIR)/image.o: image.cpp image.h graph.h kovtun.h aexpand.h TRWBP.h block.h energy.h includes.h typePotts.h MRFEnergy.h
#	$(CCOMP) $(CFLAGS) -c -o $(OBJ_DIR)/image.o image.cpp

