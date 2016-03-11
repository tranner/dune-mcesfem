set(EXTRA_DIST  parameter *.dgf)

set(DGFFILES ball.dgf circle.dgf hemi-sphere.dgf sphere.dgf spiral-2d.dgf unitcube-1d.dgf unitcube-2d.dgf ball2.dgf sphere2.dgf )

dune_symlink_to_source_files(FILES parameter ${DGFFILES})
