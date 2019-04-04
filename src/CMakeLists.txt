cmake_minimum_required(VERSION 3.0)

#set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Boost COMPONENTS program_options filesystem log iostreams REQUIRED)

set(SOURCE_FILES
        cli/command_line.cpp cli/command_line.h
        cli/exceptions.cpp cli/exceptions.h
        core/db_builder.cpp core/db_builder.h
        core/seq_traits.cpp core/seq_traits.h
        core/tree/phylo_tree.cpp core/tree/phylo_tree.h
        core/pp_matrix/proba_matrix.cpp core/pp_matrix/proba_matrix.h
        core/phylo_kmer.h
        core/pp_matrix/phyml.h core/pp_matrix/phyml.cpp
        core/kmer.h
        core/pp_matrix/branch_entry.h core/pp_matrix/branch_entry.cpp
        core/pp_matrix/row.h
        core/phylo_kmer.cpp
        core/pp_matrix/branch_entry_view.cpp core/pp_matrix/branch_entry_view.h
        utils/meta.h
        utils/file_io.cpp utils/file_io.h
        utils/algorithm.h core/phylo_kmer_db.cpp core/phylo_kmer_db.h)

include_directories(
        ./
        ../third-party/abseil-cpp
        ../third-party/fast-cpp-csv-parser
        ../third-party/range-v3/include
        #../third-party/parallel-hashmap
)

add_library(rappas-core ${SOURCE_FILES})
target_link_libraries(rappas-core
        ${Boost_LIBRARIES}
        absl::strings
        Threads::Threads
)