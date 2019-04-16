cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(builder VERSION 0.1.0 LANGUAGES CXX)

######################################################################################################
# Dependencies
find_package(Boost REQUIRED COMPONENTS program_options filesystem iostreams)
find_package(Threads REQUIRED)
find_package(rappas REQUIRED)

######################################################################################################
# Application target and properties
add_executable(builder
        # cli
        src/cli/command_line.cpp
        src/cli/exceptions.cpp

        # database builder
        src/build/db_builder.cpp

        # phylogenetic tree
        src/tree/phylo_tree.cpp

        # posterior probability matrix
        src/pp_matrix/proba_matrix.cpp
        src/pp_matrix/phyml.cpp
        src/pp_matrix/branch_entry.cpp
        src/pp_matrix/branch_entry_view.cpp

        # utils
        src/utils/file_io.cpp

        src/main.cpp
        )

target_include_directories(builder
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/
            ${CMAKE_CURRENT_SOURCE_DIR}/src/third-party/fast-cpp-csv-parser/
            ${CMAKE_CURRENT_SOURCE_DIR}/src/third-party/skarupke
)

target_link_libraries(builder
        PRIVATE
            Boost::program_options
            Boost::filesystem
            Boost::iostreams
            rappas::core
        )

# Turn on the warnings and treat them as errors
target_compile_options(builder
        PRIVATE
            -Wall -Wextra -Wpedantic -Werror
        )

target_compile_features(builder
        PUBLIC
            cxx_std_17)