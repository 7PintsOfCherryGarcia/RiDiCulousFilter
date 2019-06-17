cmake_minimum_required (VERSION 3.5)
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -g ")
project(RiDiCulousFilter)

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}"
                    ".")

# add the executable
add_executable(RiDiCulous src/main.c
                          src/main_count.c
                          src/main_gc.c
                          src/common.c)
#)
find_library(ZLIB_LIBRARY
             NAMES z)
#find_library(MLIB_LIBRARY
#             NAMES m)


get_cmake_property(_variableNames VARIABLES)

target_link_libraries(RiDiCulous
                      ${ZLIB_LIBRARY})
                      #${MLIB_LIBRARY})
