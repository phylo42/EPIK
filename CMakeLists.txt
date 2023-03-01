cmake_minimum_required(VERSION 3.10)
project(rappas2 VERSION 0.1.3 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

# default build type Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

add_subdirectory(i2l)
add_subdirectory(third-party)
add_subdirectory(rappas)