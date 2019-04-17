cmake_minimum_required(VERSION 3.13)
project(rappas_build VERSION 0.1.0 LANGUAGES CXX)

# WARNING: Some subprojects have custom libraries as dependencies:
#
#       build <- core (rappas::core)
#       build <- strasser (strasser::csv_parser)
#
# find_package(<library-name>) inside the dependent CMake files will not find required libraries unless
# you explicitly install them, for example:
#
#       cmake ..
#       make
#       make install
#
# I rely on the presence of these libraries as subprojects in this top-level directory instead, so
# no installation required.
#
# There is another workaround for the same problem
#       (see "Effective CMake" by Dan Pfeifer, https://youtu.be/bsXLMQ6WgIk)
# which is not safe though:
#       (see "Do Not Redefine CMake Commands" by Craig Scott,
#       https://crascit.com/2018/09/14/do-not-redefine-cmake-commands/)
#
# So practically we should write:
#       find_package(<library-name> REQUIRED)
# in every dependent CMake file for every dependency, but I prefer to have all of the libraries
# as subprojects and be built directly from a source tree by this CMake file.

add_subdirectory(core)
add_subdirectory(strasser)
add_subdirectory(build)