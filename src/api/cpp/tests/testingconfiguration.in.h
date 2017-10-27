//
// Created by Eric Irrgang on 10/26/17.
//

#ifndef GROMACS_TESTINGCONFIGURATION_H
#define GROMACS_TESTINGCONFIGURATION_H

namespace gmxapi
{

namespace // anonymous
{

// Todo: Need to set up a test fixture...
const std::string sample_tprfilename = "${CMAKE_CURRENT_BINARY_DIR}/topol.tpr";

} // end anonymous namespace

} // end namespace gmxapi


#endif //GROMACS_TESTINGCONFIGURATION_H
