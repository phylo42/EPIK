cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

set(ENABLE_OMP ON)
set(ENABLE_SSE OFF)
set(ENABLE_AVX2 OFF)
set(ENABLE_AVX512 OFF)

find_package(RapidJSON REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem)

if(ENABLE_OMP)
    find_package(OpenMP REQUIRED)
    message(STATUS "EPIK: OpenMP support ENABLED")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DEPIK_OMP")
else()
    message(STATUS "EPIK: OpenMP support DISABLED")
endif()

if(ENABLE_SSE)
    message(STATUS "EPIK: SSE support ENABLED")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DEPIK_SSE")
elseif(ENABLE_AVX2)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        message(STATUS "EPIK: AVX2 support ENABLED")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DEPIK_AVX2")
    endif()
elseif(ENABLE_AVX512)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        message(STATUS "EPIK: AVX2 support ENABLED")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DEPIK_AVX512")
    endif()
else()
    message(STATUS "EPIK: Vectorization DISABLED")
endif()

message(STATUS "RapidJSON: " ${RAPIDJSON_INCLUDE_DIRS})
# RapidJSON cmake scripts are different between versions
set(RapidJSON_INCLUDES ${RAPIDJSON_INCLUDE_DIRS} ${RapidJSON_INCLUDE_DIR})

set(SOURCES
        include/epik/intrinsic.h
        include/epik/jplace.h src/epik/jplace.cpp
        include/epik/place.h src/epik/place.cpp
        src/epik/main.cpp
)

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
            Boost::filesystem
            indicators::indicators
            cxxopts::cxxopts
)

if(ENABLE_OMP)
    target_link_libraries(epik-dna PRIVATE OpenMP::OpenMP_CXX)
endif()

if(ENABLE_AVX2)
    # Add compiler flags for AVX2
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(epik-dna PRIVATE -mavx2)
    endif()
endif()

if(ENABLE_AVX512)
    # Add compiler flags for AVX-512
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(epik-dna PRIVATE -mavx512f -mavx512cd)
    endif()
endif()

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
            ${RapidJSON_INCLUDES}
        )

target_link_libraries(epik-aa
        PRIVATE
        i2l::aa
        Boost::filesystem
        indicators::indicators
        cxxopts::cxxopts
        )

if(ENABLE_OMP)
    target_link_libraries(epik-aa PRIVATE OpenMP::OpenMP_CXX)
endif()


if(ENABLE_AVX512)
    # Add compiler flags for AVX-512
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(epik-aa PRIVATE -mavx512f -mavx512cd)
    endif()
endif()

target_compile_options(epik-aa
        PRIVATE
        -Wall -Wextra -Wpedantic
        )

target_compile_features(epik-aa
        PUBLIC
        cxx_std_17)

