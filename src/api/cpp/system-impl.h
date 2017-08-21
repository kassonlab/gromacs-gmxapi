#ifndef GMXAPI_SYSTEM_IMPL_H
#define GMXAPI_SYSTEM_IMPL_H

#include "gmxapi/system.h"

namespace gmxapi
{

class Runner;
class Atoms;
class MDInput;

class System::Impl
{
    public:
        Impl();
        explicit Impl(std::string tprfile);
        ~Impl();
//        int run();
//        std::unique_ptr<Atoms> atoms();
//        void setAtoms(const Atoms& atoms);
    private:
//        std::unique_ptr<Runner> runner_;
//        std::unique_ptr<Atoms>  atoms_;

};

class System::Builder
{
    public:
        Builder();
        ~Builder();

//    Builder& defaultContext();
//
//    Builder& runner(const MDInput& inputrec);
//
//    Builder& structure(const MDInput& inputrec);
//
//    Builder& topology(const MDInput& inputrec);
//
//    Builder& mdEngine(const MDInput& inputrec);

        std::unique_ptr<System> build();
    private:
        std::unique_ptr<System> system_;
};

}      // end namespace gmxapi

#endif // header guard
