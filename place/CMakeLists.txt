cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

######################################################################################################
# Application target and properties
add_executable(place_dna "")
set_target_properties(place_dna PROPERTIES OUTPUT_NAME rappas-place-dna)

target_sources(place_dna
        PRIVATE
            src/main.cpp
            src/place/place.cpp
            src/place/jplace.cpp
        )

######################################################################################################
find_package(RapidJSON REQUIRED)
target_include_directories(place_dna
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/
            ${RapidJSON_INCLUDE_DIRS}
        )

# WARNING: some dependencies are not listed here. They are supposed to be compiled
# and built together with this project by the top-level CMake file
find_package(Boost REQUIRED COMPONENTS program_options filesystem iostreams)
find_package(OpenMP REQUIRED)
find_package(QuadMath REQUIRED)

target_link_libraries(place_dna
    PRIVATE
        xpas::dna
        Boost::program_options
        Boost::filesystem
        Boost::iostreams
        OpenMP::OpenMP_CXX
        quadmath
)

# Turn on the warnings and treat them as errors
target_compile_options(place_dna
    PRIVATE
        -Wall -Wextra -Werror
)

target_compile_features(place_dna
    PUBLIC
        cxx_std_17)