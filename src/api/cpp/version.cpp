#include "gmxapi/version.h"

namespace gmxapi
{

// The values returned by these functions are compiled in when the library is built and should not
// conflict with symbols defined in a different version of the public headers that a client may
// have compiled against.
version_t Version::major()
{
    return GMXAPI_MAJOR;
}

version_t Version::minor()
{
    return GMXAPI_MINOR;
}

version_t Version::patch()
{
    return GMXAPI_PATCH;
}

std::string Version::release()
{
    return GMXAPI_RELEASE;
}

bool Version::has_feature(const std::string& featurename)
{
    // For features introduced after version 1.0, we can consult a map somewhere.
    (void)featurename;
    return false;
}

bool Version::is_at_least(version_t major, version_t minor, version_t patch)
{
    if (Version::major() >= major && Version::minor() >= minor && Version::patch() >= patch)
    {
        return true;
    }
    else
    {
        return false;
    };
}

} // end namespace gmxapi
