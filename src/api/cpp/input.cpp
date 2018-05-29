//
// Created by Eric Irrgang on 5/21/18.
//

#include <functional>
#include <map>
#include "gromacs/compat/make_unique.h"
#include "gmxapi/session/input.h"

namespace gmxapi {
namespace session {

/*!
 * \brief Variadic containe to hold maps of names to functors, with a type-based map to elements in the container.
 *
 * \tparam Ts List of element types (type keys).
 *
 * Container elements are std::map<std::string, std::function<void(T)>> and are retrieved with TypeMapMap::get(T).
 */
template<class... Ts> struct TypeMapMap{};

template<class T> struct TypeMapMap<T>
{
    std::map<std::string, std::function<void (T&)>> here;

    std::map<std::string, std::function<void (T&)>>* get(T&&)
    {
        return &here;
    };
};

template<class T, class... Ts> struct TypeMapMap<T,Ts...>
{
    std::map<std::string, std::function<void (T&)>>  here;
    TypeMapMap<Ts...> tail;

    // Something is wrong with this template or my understanding of templates.
    // U is getting deduced as a reference, but I thought the signature of get would ensure U is deduced with references removed...
    template<typename U>
    std::map<std::string, std::function<void(typename std::remove_reference<U>::type &)>>* get(U&& data)
    {
        return tail.get(std::forward<typename std::remove_reference<U>::type>(data));
    };

    std::map<std::string, std::function<void (T&)>>* get(T&&)
    {
        return &here;
    };

    std::map<std::string, std::function<void (T&)>>* get(T&)
    {
        return &here;
    };
};

class Input::Impl {
    public:
        static Input create()
        {
            std::unique_ptr<Input::Impl> implementation = gmx::compat::make_unique<Input::Impl>();
            Input output{std::move(implementation)};
            return output;
        }

        /*!
         * \brief Set the named output with the provided data.
         *
         * \tparam T Type of data output.
         * \param inputName Registered output name.
         * \param data Data of the same type as was registered for the named output.
         * \return success if the setter was called, failure if the named output is not found.
         *
         * \todo Throw an exception indicating an API bug if the type is not implemented since we are already screening
         * at the public interface.
         */
        template<typename T>
        gmxapi::Status get(const std::string &inputName,
                           T& data)
        {
            gmxapi::Status success{false}; // unsuccessful unless proven otherwise.
            try
            {
                // Not sure yet what errors occur if type of data is not in mapMap.
                std::map<std::string, std::function<void (T &)> > * getterMap = mapMap.get(data);
                auto functor = getterMap->at(inputName);
                if (functor)
                {
                    functor(data);
                    success = true;
                }
                else
                {
                    // There are no consumers registered for this output. This is not an error.
                    success = true;
                }
            }
            catch (const std::out_of_range& oor)
            {
                // The client code has attempted to access an output or type they did not register. This is a client programmer error.
            }
            return success;
        }


        // registerInput is not templated because we would want to control instantiation anyway.
        gmxapi::Status registerInput(const std::string& outputName, std::function<void (bool&)>&& functor)
        {
            // What happens if T has no default constructor?
            (*mapMap.get(bool()))[outputName] = std::move(functor);
            return gmxapi::Status(true);
        }
        gmxapi::Status registerInput(const std::string& outputName, std::function<void (double&)>&& functor)
        {
            // What happens if T has no default constructor?
            (*mapMap.get(double()))[outputName] = std::move(functor);
            return gmxapi::Status(true);
        }

        TypeMapMap<bool, double> mapMap;
};


Input::Input(std::unique_ptr<Input::Impl> &&implementation) :
    impl_{std::move(implementation)}
{
}

gmxapi::Status Input::get(const std::string &inputName,
                                 bool& data)
{
    return impl_->get(inputName, data);
}

gmxapi::Status Input::get(const std::string &inputName,
                                 double& data)
{
    return impl_->get(inputName, data);
}


Input::~Input() = default;

} // end namespace gmxapi::session
} // end namespace gmxapi
