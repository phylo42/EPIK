cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

######################################################################################################
# Dependencies
# WARNING: some dependencies are not listed here. See the top-level CMake file
# for more details
find_package(Boost REQUIRED COMPONENTS program_options filesystem iostreams)
find_package(OpenMP REQUIRED)

# The code and the libraries are the same for all targets.
# The only difference is the version of xpas we link against (xpas::dna, xpas::aa)
set(LINK_LIBRARIES Boost::program_options
        Boost::filesystem
        Boost::iostreams
        OpenMP::OpenMP_CXX
        strasser::csv_parser
        xpas::utils)

set(SOURCES src/cli/command_line.cpp
        src/cli/exceptions.cpp

        # database builder
        src/build/db_builder.cpp src/build/db_builder.h
        src/build/alignment.cpp src/build/alignment.h

        # posterior probability matrix
        src/pp_matrix/proba_matrix.cpp src/pp_matrix/proba_matrix.h
        src/pp_matrix/node_entry.cpp src/pp_matrix/node_entry.h
        src/pp_matrix/node_entry_view.cpp src/pp_matrix/node_entry_view.h
        src/pp_matrix/ar.cpp src/pp_matrix/ar.h

        src/utils/utils.h
        src/utils/io.cpp src/utils/io.h

        src/main.cpp)

######################################################################################################
# Application target and properties
add_executable(build_dna "")
set_target_properties(build_dna PROPERTIES OUTPUT_NAME rappas-build-dna)

target_sources(build_dna PRIVATE ${SOURCES})

target_include_directories(build_dna
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/
)

target_link_libraries(build_dna
        PRIVATE
            ${LINK_LIBRARIES}
            xpas::dna
        )

# Turn on the warnings and treat them as errors
target_compile_options(build_dna
    PRIVATE
        -Wall -Wextra #-Werror #-Wpedantic
    )

target_compile_features(build_dna
    PUBLIC
        cxx_std_17)

######################################################################################################
# Amino acids
add_executable(build_aa "")
set_target_properties(build_aa PROPERTIES OUTPUT_NAME rappas-build-aa)

target_sources(build_aa PRIVATE ${SOURCES})

target_include_directories(build_aa
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/
        )

target_link_libraries(build_aa
        PRIVATE
            ${LINK_LIBRARIES}
            xpas::aa
        )

# Turn on the warnings and treat them as errors
target_compile_options(build_aa
        PRIVATE
        -Wall -Wextra #-Werror #-Wpedantic
        )

target_compile_features(build_aa
        PUBLIC
        cxx_std_17)

######################################################################################################
# Amino acids, positions included

add_executable(build_aa_pos "")
set_target_properties(build_aa_pos PROPERTIES OUTPUT_NAME rappas-build-aa-pos)

target_sources(build_aa_pos PRIVATE ${SOURCES})

target_include_directories(build_aa_pos
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/
        )

target_link_libraries(build_aa_pos
        PRIVATE
        ${LINK_LIBRARIES}
        xpas::aa_pos
        )

# Turn on the warnings and treat them as errors
target_compile_options(build_aa_pos
        PRIVATE
        -Wall -Wextra #-Werror #-Wpedantic
        )

target_compile_features(build_aa_pos
        PUBLIC
        cxx_std_17)