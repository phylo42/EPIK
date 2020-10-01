cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

######################################################################################################
# Dependencies
# WARNING: some dependencies are not listed here. See the top-level CMake file
# for more details
find_package(Boost REQUIRED COMPONENTS program_options filesystem iostreams)
find_package(OpenMP REQUIRED)

######################################################################################################
# Application target and properties
add_executable(build_dna "" src/utils/utils.h)
set_target_properties(build_dna PROPERTIES OUTPUT_NAME rappas-build-dna)

target_sources(build_dna
        PRIVATE
            # cli
            src/cli/command_line.cpp
            src/cli/exceptions.cpp

            # database builder
            src/build/db_builder.cpp

            # posterior probability matrix
            src/pp_matrix/proba_matrix.cpp
            src/pp_matrix/node_entry.cpp
            src/pp_matrix/node_entry_view.cpp
        src/pp_matrix/ar.cpp

            src/main.cpp
        )

target_include_directories(build_dna
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/
)

target_link_libraries(build_dna
        PRIVATE
            Boost::program_options
            Boost::filesystem
            Boost::iostreams
            OpenMP::OpenMP_CXX
            strasser::csv_parser
            xpas::dna
            xpas::utils
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
add_executable(build_aa "" src/utils/utils.h)
set_target_properties(build_aa PROPERTIES OUTPUT_NAME rappas-build-aa)

target_sources(build_aa
        PRIVATE
        # cli
        src/cli/command_line.cpp
        src/cli/exceptions.cpp

        # database builder
        src/build/db_builder.cpp

        # posterior probability matrix
        src/pp_matrix/proba_matrix.cpp
        src/pp_matrix/node_entry.cpp
        src/pp_matrix/node_entry_view.cpp
        src/pp_matrix/ar.cpp

        src/main.cpp
        )

target_include_directories(build_aa
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/
        )

target_link_libraries(build_aa
        PRIVATE
        Boost::program_options
        Boost::filesystem
        Boost::iostreams
        OpenMP::OpenMP_CXX
        strasser::csv_parser
        xpas::aa
        xpas::utils
        )

# Turn on the warnings and treat them as errors
target_compile_options(build_aa
        PRIVATE
        -Wall -Wextra #-Werror #-Wpedantic
        )

target_compile_features(build_aa
        PUBLIC
        cxx_std_17)