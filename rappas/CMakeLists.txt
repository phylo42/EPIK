cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

set(SOURCES
        src/main.cpp
        src/rappas/place.cpp src/rappas/place.h
        src/rappas/jplace.cpp src/rappas/jplace.h)

######################################################################################################
# Application target and properties
add_executable(rappas2-dna "")

target_sources(rappas2-dna
        PRIVATE
            ${SOURCES})

find_package(RapidJSON REQUIRED)
target_include_directories(rappas2-dna
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
            xcl::dna
            Boost::program_options
            Boost::filesystem
            Boost::iostreams
            OpenMP::OpenMP_CXX
            quadmath
)

# Turn on the warnings and treat them as errors
target_compile_options(rappas2-dna
        PRIVATE
            -Wall -Wextra -Werror
)


target_compile_features(rappas2-dna
        PUBLIC
            cxx_std_17)

######################################################################################################

#add_executable(rappas2-aa "")

#target_sources(rappas2-aa
#        PRIVATE
#        ${SOURCES})

#target_include_directories(rappas2-aa
#        PRIVATE
#            ${CMAKE_CURRENT_SOURCE_DIR}/src/
#            ${RapidJSON_INCLUDE_DIRS}
#        )

#target_link_libraries(rappas2-aa
#        PRIVATE
#        xcl::aa
#        Boost::program_options
#        Boost::filesystem
#        Boost::iostreams
#        OpenMP::OpenMP_CXX
#        quadmath
#        )

#target_compile_options(rappas2-aa
#        PRIVATE
#        -Wall -Wextra -Werror
#        )

#target_compile_features(rappas2-aa
#        PUBLIC
#        cxx_std_17)
#
