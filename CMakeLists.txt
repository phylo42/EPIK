cmake_minimum_required(VERSION 3.10)
project(rappas_build VERSION 0.1.0 LANGUAGES CXX)

# default build type Release
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# define a hash map
if(NOT HASH_MAP)
    set(HASH_MAP "USE_SKA_BYTELL_HASH_MAP")
endif()

add_subdirectory(core)
add_subdirectory(strasser)
add_subdirectory(build)