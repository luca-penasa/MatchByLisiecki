cmake_minimum_required(VERSION 2.8)
project(Match)


file( GLOB_RECURSE header_list *.hh)
file( GLOB_RECURSE source_list *.cc)

add_executable(match ${header_list} ${source_list})
