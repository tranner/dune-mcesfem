# include not needed for CMake
# include $(top_srcdir)/am/global-rules
set(CLEANFILES  *.log *.vtu *.pvd *.pvtu *.series *.txt)

dune_symlink_to_source_files(FILES errorplot.py)
