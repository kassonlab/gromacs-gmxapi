//
// Created by Eric Irrgang on 5/21/18.
//

#include "gmxapi/context/outputstream.h"

namespace gmxapi {
namespace context {

class OutputStream::Impl {
};

void OutputStream::set(std::string outputName,
                       bool data)
{}

OutputStream::~OutputStream() = default;

} // end namespace gmxapi::context
} // end namespace gmxapi
