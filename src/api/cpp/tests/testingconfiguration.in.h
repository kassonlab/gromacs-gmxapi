//
// Created by Eric Irrgang on 10/26/17.
//

#ifndef GROMACS_TESTINGCONFIGURATION_H
#define GROMACS_TESTINGCONFIGURATION_H

#include <string>
#include <vector>

namespace gmxapi
{

namespace testing
{

// Todo: Need to set up a test fixture...
static const std::string sample_tprfilename = "${CMAKE_CURRENT_BINARY_DIR}/topol.tpr";

static const std::vector<std::string> mdArgs{"-ntomp", "1"};

} // end namespace gmxapi::testing

} // end namespace gmxapi


#endif //GROMACS_TESTINGCONFIGURATION_H
