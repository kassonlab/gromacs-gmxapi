//
// Created by Eric Irrgang on 5/21/18.
//

#ifndef GROMACS_OUTPUTSTREAM_H
#define GROMACS_OUTPUTSTREAM_H

#include <memory>
#include <string>

namespace gmxapi {
namespace context {

class OutputStream
{
    public:
        ~OutputStream();

        /*!
         * \brief Set data for a registered output stream.
         *
         * \param outputName Registered name of the output port
         * \param data data to set with the registered output handler.
         *
         * We should not use a template here to handle the different data types because the template might be expressed
         * with different munged symbol names by different compilers. But we want this interface to be extensible, so
         * we need to consider how to allow flexible types. We could wrap all data in a gmxapi::Data wrapper or something,
         * but that makes calling the set() method more cumbersome in the client code.
         *
         * What we could do, though, is to provide a template function as a helper that is compiled in the client code
         * and just makes it easier to make the correct call here. Then a gmxapi::Data wrapper wouldn't be cumbersome.
         */
        void set(std::string outputName, bool data);
        //void set(std::string outputName, someOtherType

        void set(const char* outputName, bool data)
        {
            this->set(std::string(outputName), data);
        }

    private:
        // Private implementation class
        class Impl;
        // opaque pointer to implementation.
        std::unique_ptr<Impl> impl_;
};

} // end namespace gmxapi::context
} //end namespace gmxapi


#endif //GROMACS_OUTPUTSTREAM_H
