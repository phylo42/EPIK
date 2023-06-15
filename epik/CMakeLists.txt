cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

find_package(RapidJSON REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options filesystem iostreams)
find_package(OpenMP REQUIRED)
find_package(QuadMath REQUIRED)


set(SOURCES
        src/epik/main.cpp
        include/epik/place.h src/epik/place.cpp
        include/epik/jplace.h src/epik/jplace.cpp)

set(LIBRARIES
        Boost::program_options
        Boost::filesystem
        Boost::iostreams
        OpenMP::OpenMP_CXX
        quadmath)

######################################################################################################
# Application target and properties
add_executable(epik-dna "")
target_sources(epik-dna PRIVATE ${SOURCES})

target_include_directories(epik-dna
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/
        ${CMAKE_CURRENT_SOURCE_DIR}/include/
        ${RapidJSON_INCLUDE_DIR}
        )

target_link_libraries(epik-dna
        PRIVATE
            i2l::dna
            ${LIBRARIES}
        )

# Turn on the warnings and treat them as errors
target_compile_options(epik-dna
        PRIVATE
            -Wall -Wextra -Wpedantic
)


target_compile_features(epik-dna
        PUBLIC
            cxx_std_17)

######################################################################################################

add_executable(epik-aa "")

target_sources(epik-aa
        PRIVATE
        ${SOURCES})

target_include_directories(epik-aa
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/
            ${CMAKE_CURRENT_SOURCE_DIR}/include/
            ${RapidJSON_INCLUDE_DIR}
        )

target_link_libraries(epik-aa
        PRIVATE
            i2l::aa
            ${LIBRARIES}
        )

target_compile_options(epik-aa
        PRIVATE
        -Wall -Wextra -Wpedantic
        )

target_compile_features(epik-aa
        PUBLIC
        cxx_std_17)

