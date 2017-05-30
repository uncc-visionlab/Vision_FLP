cmake_minimum_required(VERSION 2.8)

#SET(CMAKE_BUILD_TYPE distribution)
#SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")
#SET(CMAKE_C_FLAGS_DISTRIBUTION "-O3")

#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/lib)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/bin)


project(MRFCode)
subdirs(src/bp-vision)
subdirs(src/MRF2.2)
#subdirs(src/TRW_S-v1.3)
subdirs(src/eff_mlabel_ver1)

configure_file(src/bp-vision/penguin.pgm 
    ${CMAKE_BINARY_DIR}/build/bin/penguin.pgm COPYONLY)
configure_file(src/bp-vision/tsukuba1.pgm 
    ${CMAKE_BINARY_DIR}/build/bin/tsukuba1.pgm COPYONLY)
configure_file(src/bp-vision/tsukuba2.pgm
    ${CMAKE_BINARY_DIR}/build/bin/tsukuba2.pgm COPYONLY)
configure_file(src/bp-vision/dimage.dat
    ${CMAKE_BINARY_DIR}/build/bin/images/dimage.dat COPYONLY)

configure_file(src/eff_mlabel_ver1/images/datacost.txt 
    ${CMAKE_BINARY_DIR}/build/bin/images/datacost.txt COPYONLY)
configure_file(src/eff_mlabel_ver1/images/tsukuba_l.ppm 
    ${CMAKE_BINARY_DIR}/build/bin/images/tsukuba_l.ppm COPYONLY)
configure_file(src/eff_mlabel_ver1/images/tsukuba_r.ppm
    ${CMAKE_BINARY_DIR}/build/bin/images/tsukuba_r.ppm COPYONLY)

