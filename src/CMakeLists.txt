cmake_minimum_required(VERSION 3.0)
project(rappas-cpp)

find_package(Boost COMPONENTS program_options filesystem log iostreams REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(SOURCE_FILES
        main.cpp
        core/db_builder.cpp core/db_builder.h
        cli/command_line.cpp cli/command_line.h
        cli/exceptions.cpp cli/exceptions.h
        ar/evo_model.cpp ar/evo_model.h
        ar/bindep_wrapper.cpp ar/bindep_wrapper.h
        core/alignment.cpp core/alignment.h
        core/fasta.cpp core/fasta.h
        core/seq_traits.cpp core/seq_traits.h core/newick.cpp core/newick.h)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
add_executable(rappas-cpp ${SOURCE_FILES})
target_link_libraries(rappas-cpp ${Boost_LIBRARIES})
