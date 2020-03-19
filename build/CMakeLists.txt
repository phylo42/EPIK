cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

######################################################################################################
# Dependencies
# WARNING: some dependencies are not listed here. See the top-level CMake file
# for more details
find_package(Boost REQUIRED COMPONENTS program_options filesystem iostreams)
find_package(OpenMP REQUIRED)

######################################################################################################
# Application target and properties
add_executable(build_n "")
set_target_properties(build_n PROPERTIES OUTPUT_NAME rappas-buildn)

target_sources(build_n
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
            src/pp_matrix/phyml.cpp

            src/main.cpp
        )

target_include_directories(build_n
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/
)

target_link_libraries(build_n
        PRIVATE
            Boost::program_options
            Boost::filesystem
            Boost::iostreams
            OpenMP::OpenMP_CXX
            strasser::csv_parser
            rappas::core_n
            rappas::utils
        )

# Turn on the warnings and treat them as errors
target_compile_options(build_n
        PRIVATE
            -Wall -Wextra #-Werror #-Wpedantic
        )

target_compile_features(build_n
        PUBLIC
            cxx_std_17)