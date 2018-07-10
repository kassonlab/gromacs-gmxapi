# This list file provides the gmxapi::gmxapi cmake module. gmxapi::tests is the lone
# submodule, but subdirectories include the public headers.

# Define the gmxapi library target. Dependent targets can use properties of the target
# for include directories and header files.

# Configure / copy header files from the source to the build directory for clarity and packaging
# high-level client headers
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gmxapi)
configure_file(gmxapi/context.h gmxapi)
configure_file(gmxapi/exceptions.h gmxapi)
configure_file(gmxapi/gromacsfwd.h gmxapi)
configure_file(gmxapi/gmxapi.h gmxapi)
configure_file(gmxapi/md.h gmxapi)
#configure_file(gmxapi/runner.h gmxapi)
configure_file(gmxapi/session.h gmxapi)
configure_file(gmxapi/status.h gmxapi)
configure_file(gmxapi/system.h gmxapi)
configure_file(gmxapi/version.in.h gmxapi/version.h)

# headers for API development and extension

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gmxapi/md)
#configure_file(gmxapi/md/runnerstate.h gmxapi/md)
configure_file(gmxapi/md/mdmodule.h gmxapi/md)
configure_file(gmxapi/md/mdsignals.h gmxapi/md)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gmxapi/session)
configure_file(gmxapi/session/outputstream.h gmxapi/session)
configure_file(gmxapi/session/input.h gmxapi/session)

# Add to install target: copy the public API headers from the source directory
# into the gmxapi header subdirectory.
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gmxapi
        DESTINATION include
        FILES_MATCHING PATTERN "*.h")

# Using $<TARGET_OBJECTS:mdrun_objlib> doesn't work in all environments when linking a library.
# The object files in programs/mdrun/ are compiled without -fPIC by default and so may not be
# linkable into a shared object. I had to tell cmake to build mdrun_objlib relocatable...
add_library(gmxapi SHARED $<TARGET_OBJECTS:mdrun_objlib>)

target_include_directories(gmxapi PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include>
    )

target_include_directories(gmxapi PRIVATE
                           ${CMAKE_CURRENT_BINARY_DIR}
                           ${CMAKE_CURRENT_SOURCE_DIR}
                           )

# libgmxapi source files. The public headers are picked up automatically by CMake
# for the build system dependencies, but since they aren't listed as sources, older Cmake
# might not get triggered by updates to the public headers. This might be reason to
# use an INTERFACE library target, but since we call configure_file() above to get everything
# in one place for packaging, CMake is already making sure the public headers are
# up-to-date. However, an update to an existing public header that adds a dependency
# on another header would _not_ be expected to trigger a rebuild of the gmxapi target
# in the current setup. This is an argument against headers including headers, since
# to create a target (including an INTERFACE target) that has all of the headers
# would trigger the build of any dependent target even if it does not depend on a
# given header in the public interface. This isn't a huge deal, and a public header
# change probably implies an important API change, but it might also just be a result
# of a documentation change or something that does not affect the part of the ABI
# used by a given client. If it were easy to use the sources from an interface target
# to install the public headers, I would think we should use it, but it is not, so
# I won't.
# Update: more recent versions of CMake seem to do more rigorous dependency chasing,
# but listing all source files helps some code introspection tools work better.
target_sources(
        gmxapi PRIVATE
        context.cpp
        gmxapi.cpp
        input.cpp
        md-impl.h
        md.cpp
        mdmodule.cpp
        mdsignals.cpp
        outputstream.cpp
        session-impl.h
        session.cpp
        status.cpp
        system.cpp
        system-impl.h
        version.cpp
        workflow.h
        workflow-impl.h
        workflow.cpp
        )

# Ref. https://cmake.org/Wiki/CMake_RPATH_handling
# use, i.e. don't skip the full RPATH for the build tree
set_target_properties(gmxapi PROPERTIES SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set_target_properties(gmxapi PROPERTIES BUILD_WITH_INSTALL_RPATH FALSE)

# duplicates link path when installed in LIB_INSTALL_DIR with libgromacs
# but without it, targets referencing libgromacs symbols when importing only this target will
# not find libgromacs the way OS X has set the RPATH to @executable_path=../lib
# To set the @rpath to resolve to the directory in which the libgmxapi shared
# object is living, we would set it to @loader_path, but I don't know of a CMake
# directive that produces that effect. We may just need a platform check.
# I suppose the question is what/how/why does CMake set rpath to executable_path instead
# of rpath...
set_target_properties(gmxapi PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}")

set_target_properties(gmxapi PROPERTIES
                      INSTALL_RPATH_USE_LINK_PATH TRUE
                      SOVERSION ${GMXAPI_MAJOR}
                      VERSION ${GMXAPI_RELEASE}
                      )

target_link_libraries(gmxapi PRIVATE libgromacs)

# Install the gmxapi target and simultaneously define the export target for
# which CMake will create a helper file. Specify the directory for clients to
# add to their include path to be able to `#include "gmxapi/some_header.h"`
install(TARGETS gmxapi
        EXPORT gmxapi
        LIBRARY DESTINATION ${LIB_INSTALL_DIR}
        ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
        RUNTIME DESTINATION ${LIB_INSTALL_DIR}
        INCLUDES DESTINATION include
        )


# Create the CMake exports file to help other projects build against libgmxapi
# as a CMake import target Gromacs::gmxapi.
install(EXPORT gmxapi
        NAMESPACE Gromacs::
        DESTINATION share/cmake/gmxapi/
        )
add_library(Gromacs::gmxapi ALIAS gmxapi )

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/gmxapi-config-version.cmake"
    VERSION ${GMXAPI_RELEASE}
    COMPATIBILITY ExactVersion
    # Every dev release will break things for a while. We can assert backwards compatibility closer to 1.0
    # COMPATIBILITY SameMajorVersion
)
configure_package_config_file(
    cmake/gmxapi-config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/gmxapi-config.cmake"
    INSTALL_DESTINATION share/cmake/gmxapi/
    PATH_VARS LIB_INSTALL_DIR
)
install(
    FILES
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/gmxapi-config-version.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/gmxapi-config.cmake
    DESTINATION share/cmake/gmxapi/
)


if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
