cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(GITHOME_MSU $ENV{GITHOME_MSU})

project(EoSTest CXX)

add_custom_target(extern_msu_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)
add_custom_target(extern_msu_sampler COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/msu_sampler/software)
add_custom_target(extern_msu_eos COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/msu_eos/software)

set(CMAKE_CXX_FLAGS "-O2 -Wall -std=c++17 -mfpmath=sse")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/msu_sampler/software/include ${GITHOME_MSU}/msu_eos/software/include)

link_directories(${GITHOME_MSU}/msu_sampler/software/lib ${GITHOME_MSU}/commonutils/software/lib ${GITHOME_MSU}/msu_eos/software/lib)

find_library(${GITHOME_MSU}/msu_sampler/software/lib msu_sampler)
find_library(${GITHOME_MSU}/commonutils/software/lib msu_commonutils)
find_library(${GITHOME_MSU}/msu_eos/software/lib msu_eos)

#######
set(SOURCE_FILES ./eostest.cc)
SET_SOURCE_FILES_PROPERTIES(eostest.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/msu_sampler/software/lib/libmsu_sampler.a)
SET_SOURCE_FILES_PROPERTIES(eostest.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/msu_eos/software/lib/libmsu_eos.a)
SET_SOURCE_FILES_PROPERTIES(eostest.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a)
	  
add_executable(eostest ${SOURCE_FILES})
add_dependencies(eostest extern_msu_sampler)
add_dependencies(eostest extern_msu_eos)
add_dependencies(eostest extern_msu_commonutils)
target_link_libraries(eostest msu_sampler msu_commonutils msu_eos ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY} ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GITHOME_MSU}/msu_eos/software/lib/libmsu_eos.a ${GITHOME_MSU}/msu_sampler/software/lib/libmsu_sampler.a)
########

#######
set(SOURCE_FILES ./eostest_vsepsilon.cc)
SET_SOURCE_FILES_PROPERTIES(eostest_vsepsilon.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/msu_sampler/software/lib/libmsu_sampler.a)
SET_SOURCE_FILES_PROPERTIES(eostest_vsepsilon.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/msu_eos/software/lib/libmsu_eos.a)
SET_SOURCE_FILES_PROPERTIES(eostest_vsepsilon.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a)
	  
add_executable(eostest_vsepsilon ${SOURCE_FILES})
add_dependencies(eostest_vsepsilon extern_msu_sampler)
add_dependencies(eostest_vsepsilon extern_msu_eos)
add_dependencies(eostest_vsepsilon extern_msu_commonutils)
target_link_libraries(eostest_vsepsilon msu_sampler msu_commonutils msu_eos ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY} ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GITHOME_MSU}/msu_eos/software/lib/libmsu_eos.a ${GITHOME_MSU}/msu_sampler/software/lib/libmsu_sampler.a)
########
	