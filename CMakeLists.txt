cmake_minimum_required(VERSION 3.10)
project(rappas_build VERSION 0.1.3 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

# default build type Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# xpas can be compiled to either store phylo k-mer positions or not.
# Comment it out if you do not want to store them.
# WARNING! Setting it to 0, false or whatever else would not help,
# xpas checks if it is DEFINED no matter what is the value.
#set(KEEP_POSITIONS 1)

# define a hash map
if(NOT DEFINED HASH_MAP)
    set(HASH_MAP "USE_TSL_ROBIN_MAP")
endif()

add_subdirectory(xpas)
add_subdirectory(strasser)
add_subdirectory(build)
add_subdirectory(place)