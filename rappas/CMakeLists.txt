cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

######################################################################################################
# Application target and properties
add_executable(rappas2-dna "")
add_executable(rappas2-aa "")


set(SOURCES
        src/main.cpp
        src/rappas/place.cpp src/rappas/place.h
        src/rappas/jplace.cpp src/rappas/jplace.h)

target_sources(rappas2-dna
        PRIVATE
            ${SOURCES})

target_sources(rappas2-aa
        PRIVATE
        ${SOURCES})

######################################################################################################
find_package(RapidJSON REQUIRED)
target_include_directories(rappas2-dna
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/
            ${RapidJSON_INCLUDE_DIRS}
        )
target_include_directories(rappas2-aa
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/
            ${RapidJSON_INCLUDE_DIRS}
        )

# WARNING: some dependencies are not listed here. They are supposed to be compiled
# and built together with this project by the top-level CMake file
find_package(Boost REQUIRED COMPONENTS program_options filesystem iostreams)
find_package(OpenMP REQUIRED)
find_package(QuadMath REQUIRED)

target_link_libraries(rappas2-dna
        PRIVATE
            xpas::dna
            Boost::program_options
            Boost::filesystem
            Boost::iostreams
            OpenMP::OpenMP_CXX
            #quadmath
)

target_link_libraries(rappas2-aa
        PRIVATE
            xpas::aa
            Boost::program_options
            Boost::filesystem
            Boost::iostreams
            OpenMP::OpenMP_CXX
            #quadmath
        )

# Turn on the warnings and treat them as errors
target_compile_options(rappas2-dna
        PRIVATE
            -Wall -Wextra #-Werror
)

set_property(TARGET rappas2-aa PROPERTY CXX_STANDARD 17)
target_compile_options(rappas2-aa
        PRIVATE
            -Wall -Wextra #-Werror
        )

set_property(TARGET rappas2-dna PROPERTY CXX_STANDARD 17)
target_compile_features(rappas2-dna
        PUBLIC
            cxx_std_17)

target_compile_features(rappas2-aa
        PUBLIC
            cxx_std_17)

