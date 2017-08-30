#include "io.h"

#include <memory>

#include "gmxapi/md.h"

namespace gmxapi
{
namespace io
{

TprFile::TprFile(const std::string &filename, const file_mode mode)
{
    if (mode != file_mode::READ)
    {
        // TODO: helpful exception and message
        throw(gmxapi::Exception());
    }
    // TODO: check file type.

    // /* Read (nearly) all data required for the simulation */
    // read_tpx_state(ftp2fn(efTPR, nfile, fnm), inputrec, state, mtop);

    // TODO: confirm file type and file accessibility.

    t_state    state;    // Integrator state
    t_inputrec inputrec; // simulation parameters
    gmx_mtop_t mtop;     // forcefield params, atom typing, molecules

    // Opens and closes TPR file.
    read_tpx_state(filename.cstr(), &inputrec, &state, &mtop);
}

std::unique_ptr<MDInput> TprFile::mdInput()
{}

TprFile::~TprFile()
{};

} // end namespace io
} // end namespace gmxapi
