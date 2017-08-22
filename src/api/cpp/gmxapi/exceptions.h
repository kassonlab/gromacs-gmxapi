#ifndef GMXAPI_EXCEPTIONS_H
#define GMXAPI_EXCEPTIONS_H
/*! \file
 * \brief Declare exception classes for external API.
 *
 * \ingroup gmxapi
 */

#include <exception>

namespace gmxapi
{

/*! \brief Base exception for gmxapi library.
 *
 * Exceptions thrown in the gmxapi namespace are descended from gmxapi::Exception
 * or there is a bug.
 *
 * \ingroup gmxapi
 */
class Exception : public std::exception
{
    public:
        const char* what() const noexcept override
        {
            return "Gromacs API error";
        };
};

}      // end namespace gmxapi

#endif // header guard
