#ifndef GMXAPI_IO_H
#define GMXAPI_IO_H
//! \cond
/*! \file
 * \brief Declare I/O functionality for API implementation helpers.
 *
 */

#include <string>

namespace gmxapi
{
namespace io
{

/// Object-oriented access to TPR files.
/*! TprFile objects are explicitly file-backed containers for simulation parameters,
 * structure and topology data. The file is opened at object creation and
 * guaranteed to be closed after object destruction.
 */
class TprFile
{
    public:
        // Should this be a more general data access mode?
        // If so, it may be more flexible to use bit masks, since write-only data makes sense in some contexts.
        // Streams are
        // append-only or read-once & non-seekable.
        enum class file_mode
        {
            READ,
            WRITE,
            APPEND
        };


        explicit TprFile(const std::string &filename, const file_mode mode);

        ~TprFile();
    private:
        const std::string filename_;
};

}      // end namespce gmxapi::io
}      // end namespce gmxapi
//! \endcond
#endif // header guard
