cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

######################################################################################################
# Application target and properties
add_executable(place_n "")
set_target_properties(place_n PROPERTIES OUTPUT_NAME rappas-placen)

target_sources(place_n
        PRIVATE
            src/main.cpp
            src/place/place.cpp
            src/place/jplace.cpp
        )

######################################################################################################
find_package(RapidJSON REQUIRED)
target_include_directories(place_n
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/
            ${RapidJSON_INCLUDE_DIRS}
        )

# WARNING: some dependencies are not listed here. They are supposed to be compiled
# and built together with this project by the top-level CMake file
find_package(Boost REQUIRED COMPONENTS program_options filesystem iostreams)
find_package(OpenMP REQUIRED)
find_package(QuadMath REQUIRED)

target_link_libraries(place_n
    PRIVATE
        xpas::dna
        Boost::program_options
        Boost::filesystem
        Boost::iostreams
        OpenMP::OpenMP_CXX
        quadmath
)

# Turn on the warnings and treat them as errors
target_compile_options(place_n
    PRIVATE
        -Wall -Wextra -Werror
)

target_compile_features(place_n
    PUBLIC
        cxx_std_17)