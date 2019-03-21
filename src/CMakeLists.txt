cmake_minimum_required(VERSION 3.0)

#set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Boost COMPONENTS program_options filesystem log iostreams REQUIRED)

set(SOURCE_FILES
        #main.cpp
        cli/command_line.cpp cli/command_line.h
        cli/exceptions.cpp cli/exceptions.h
        core/db_builder.cpp core/db_builder.h
        core/seq_traits.cpp core/seq_traits.h
        core/phylo_tree.cpp core/phylo_tree.h
        core/proba_matrix.cpp core/proba_matrix.h
        core/phylo_kmer_explorer.cpp
        core/phylo_kmer_explorer.h
        core/phylo_kmer.h
        utils/file_io.cpp utils/file_io.h core/phyml.h core/phyml.cpp utils/algorithm.h core/kmer.h)


include_directories(
        ../third-party/abseil-cpp
        ../third-party/fast-cpp-csv-parser
)

add_library(rappas-core ${SOURCE_FILES})
target_link_libraries(rappas-core
        ${Boost_LIBRARIES}
        absl::strings
        Threads::Threads
)