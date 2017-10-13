include_directories(${CMAKE_SOURCE_DIR}/src/external/gmock-1.7.0/gtest/include
                    ${GMOCK_INCLUDE_DIRS})

add_executable(harmonic-test test_harmonic.cpp)
add_test(HarmonicPotentialPlugin harmonic-test)

target_link_libraries(harmonic-test harmonicpotential libgromacs Gromacs::gmxapi
                      gmock
                      ${GMOCK_LIBRARIES})

set_tests_properties(HarmonicPotentialPlugin PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                     COMPILE_DEFINITIONS "${GMOCK_COMPILE_DEFINITIONS}"
                     COMPILE_FLAGS "${GMOCK_COMPILE_FLAGS}")
