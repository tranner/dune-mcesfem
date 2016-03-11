set(GRIDTYPE  ALUGRID_SIMPLEX)
set(GRIDDIM  3)
set(WORLDDIM  3)
set(POLORDER  2)
set(WANT_ISTL  1)

set(CLEANFILES  *.log)

set(check_PROGRAMS  main )

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DPOLORDER=${POLORDER} -DWANT_ISTL=${WANT_ISTL} -DGRIDDIM=${GRIDDIM} -DWORLDDIM=${WORLDDIM}")

add_executable(main main.cc)
add_dune_alberta_flags(main GRIDDIM ${GRIDDIM} WORLDDIM ${WORLDDIM})
target_link_dune_default_libraries( main )

set(EXTRA_DIST   main.cc deformation.hh dirichletconstraints.hh elliptic.hh femscheme.hh heat.hh heatmodel.hh heatscheme.hh model.hh probleminterface.hh rhs.hh temporalprobleminterface.hh heat2.dox)

#add_directory_test_target(_test_target)
#add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
