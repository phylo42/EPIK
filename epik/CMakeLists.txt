cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

set(ENABLE_OMP ON)

find_package(RapidJSON REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options filesystem iostreams)
find_package(QuadMath REQUIRED)
if(ENABLE_OMP)
    find_package(OpenMP REQUIRED)
    message(STATUS "OpenMP: ENABLED")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DEPIK_OMP")
else()
    message(STATUS "OpenMP: DISABLED")
endif()

message(STATUS "RapidJSON: " ${RAPIDJSON_INCLUDE_DIRS})
# RapidJSON cmake scripts are different between versions
set(RapidJSON_INCLUDES ${RAPIDJSON_INCLUDE_DIRS} ${RapidJSON_INCLUDE_DIR})

set(SOURCES
        src/epik/main.cpp
        include/epik/place.h src/epik/place.cpp
        include/epik/jplace.h src/epik/jplace.cpp)

######################################################################################################
# Application target and properties
add_executable(epik-dna "")
target_sources(epik-dna PRIVATE ${SOURCES})

target_include_directories(epik-dna
        PRIVATE
        ${RAPIDJSON_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/
        ${CMAKE_CURRENT_SOURCE_DIR}/include/
        )

target_link_libraries(epik-dna
        PRIVATE
            i2l::dna
            Boost::program_options
            Boost::filesystem
            Boost::iostreams
            quadmath
)

if(ENABLE_OMP)
    target_link_libraries(epik-dna PRIVATE OpenMP::OpenMP_CXX)
endif()

# Turn on the warnings and treat them as errors
target_compile_options(epik-dna
        PRIVATE
            -Wall -Wextra -Werror
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
            ${RapidJSON_INCLUDES}
        )

target_link_libraries(epik-aa
        PRIVATE
        i2l::aa
        Boost::program_options
        Boost::filesystem
        Boost::iostreams
        quadmath
        )
if(ENABLE_OMP)
    target_link_libraries(epik-aa PRIVATE OpenMP::OpenMP_CXX)
endif()

target_compile_options(epik-aa
        PRIVATE
        -Wall -Wextra -Werror
        )

target_compile_features(epik-aa
        PUBLIC
        cxx_std_17)

