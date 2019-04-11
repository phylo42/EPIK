cmake_minimum_required(VERSION 3.0)

set(PROJECT_NAME rappas-build)
project(${PROJECT_NAME})

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Boost COMPONENTS program_options filesystem iostreams REQUIRED)
find_package(Threads REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror")

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG)

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_DEBUG "-g")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)


set(SOURCE_FILES
        # cli
        src/cli/command_line.cpp
        src/cli/command_line.h
        src/cli/exceptions.cpp
        src/cli/exceptions.h

        # database builder
        src/build/db_builder.cpp
        src/build/db_builder.h

        # phylogenetic tree
        src/tree/phylo_tree.cpp
        src/tree/phylo_tree.h

        # posterior probability matrix
        src/pp_matrix/proba_matrix.cpp
        src/pp_matrix/proba_matrix.h
        src/pp_matrix/phyml.h
        src/pp_matrix/phyml.cpp
        src/pp_matrix/branch_entry.h
        src/pp_matrix/branch_entry.cpp
        src/pp_matrix/row.h
        src/pp_matrix/branch_entry_view.cpp
        src/pp_matrix/branch_entry_view.h

        # utils
        src/utils/meta.h
        src/utils/file_io.cpp
        src/utils/file_io.h
        src/utils/algorithm.h

        src/main.cpp
        )

include_directories(src)
include_directories(third-party/skarupke)
include_directories(third-party/fast-cpp-csv-parser)
include_directories(dependencies/rappas-core/include)

add_subdirectory(dependencies/rappas-core)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME}
        rappas-core
        ${Boost_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT})