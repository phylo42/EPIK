cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

######################################################################################################
# Application target and properties
add_executable(place_dna "")
add_executable(place_aa "")

set_target_properties(place_dna PROPERTIES OUTPUT_NAME rappas-place-dna)
set_target_properties(place_aa PROPERTIES OUTPUT_NAME rappas-place-aa)

set(SOURCES
        src/main.cpp
        src/place/place.cpp
        src/place/jplace.cpp)

target_sources(place_dna
        PRIVATE
            ${SOURCES})

target_sources(place_aa
        PRIVATE
        ${SOURCES})

######################################################################################################
find_package(RapidJSON REQUIRED)
target_include_directories(place_dna
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/
            ${RapidJSON_INCLUDE_DIRS}
        )
target_include_directories(place_aa
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

target_link_libraries(place_aa
        PRIVATE
            xpas::aa
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
target_compile_options(place_aa
        PRIVATE
            -Wall -Wextra -Werror
        )

target_compile_features(place_dna
        PUBLIC
            cxx_std_17)

target_compile_features(place_aa
        PUBLIC
            cxx_std_17)

