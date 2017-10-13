//
// Created by Eric Irrgang on 10/13/17.
//

#include "harmonicpotential.h"

namespace plugin
{

Harmonic::~Harmonic() = default;

gmx::vec3<real> Harmonic::calculateForce(gmx::vec3<real> r1,
                                         gmx::vec3<real> r2)
{
    return gmx::vec3<real>(1.f, 1.f, 1.f);
}
} // end namespace plugin
