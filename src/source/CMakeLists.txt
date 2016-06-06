set(GRIDTYPE  ALBERTAGRID)
set(GRIDDIM  2)
set(WORLDDIM  3)
set(POLORDER  1)
set(WANT_ISTL  1)

set(CLEANFILES  *.log)

set(check_PROGRAMS  main main-2d )

# SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DPOLORDER=${POLORDER} -DWANT_ISTL=${WANT_ISTL} -DGRIDDIM=${GRIDDIM} -DWORLDDIM=${WORLDDIM}")

add_executable(main main.cc)
add_dune_alberta_flags(main GRIDDIM 2 WORLDDIM 3)
target_link_dune_default_libraries( main )
target_compile_definitions( main PRIVATE ${GRIDTYPE};GRIDDIM=2;WORLDDIM=3;WANT_ISTL=${WANT_ISTL};POLORDER=${POLORDER} )

set(EXTRA_DIST   main.cc deformation.hh dirichletconstraints.hh elliptic.hh femscheme.hh heat.hh heatmodel.hh heatscheme.hh model.hh probleminterface.hh rhs.hh temporalprobleminterface.hh heat2.dox)

add_executable(main-2d main.cc)
add_dune_alberta_flags(main-2d GRIDDIM 1 WORLDDIM 2 )
target_link_dune_default_libraries( main-2d )
target_compile_definitions( main-2d PRIVATE GRIDDIM=1;WORLDDIM=2;WANT_ISTL=${WANT_ISTL};POLORDER=${POLORDER} )


#add_directory_test_target(_test_target)
#add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})

dune_symlink_to_source_files(FILES run.py 2d-run.sh)
