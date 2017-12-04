#
#add_executable(restraintpotential-test
#               manager.cpp
#               $<TARGET_OBJECTS:restraint>)

gmx_add_unit_test(RestraintTests restraintpotential-test
                  manager.cpp
                  $<TARGET_OBJECTS:restraint>)
