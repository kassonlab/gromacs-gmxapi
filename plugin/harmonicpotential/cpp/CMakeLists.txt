add_library(harmonicpotential SHARED
            harmonicpotential.h
            harmonicpotential.cpp)

target_include_directories(harmonicpotential PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:include>
                           )
target_link_libraries(harmonicpotential libgromacs)
