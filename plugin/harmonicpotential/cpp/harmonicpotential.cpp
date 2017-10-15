//
// Created by Eric Irrgang on 10/13/17.
//

#include <cmath>
#include "harmonicpotential.h"

namespace plugin
{

Harmonic::~Harmonic() = default;

gmx::vec3<real> Harmonic::calculateForce(gmx::vec3<real> r1,
                                         gmx::vec3<real> r2)
{
    // set equilibrium separation distance
    // TODO: be clearer about units
    real x0{1.0};
    // set spring constant
    // TODO: be clearer about units
    real k{1.0};
    auto r12 = r2 - r1;
    // TODO: find appropriate math header and namespace
    auto r = sqrt(dot(r12, r12));
    return gmx::vec3<real>(1.f, 1.f, 1.f);
}
} // end namespace plugin
