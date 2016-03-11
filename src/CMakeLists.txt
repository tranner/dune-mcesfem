add_subdirectory("data")
add_subdirectory("output")
add_subdirectory("source")

add_executable("dune-mcesfem" dune-mcesfem.cc)
target_link_dune_default_libraries("dune-mcesfem")
