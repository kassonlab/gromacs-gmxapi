#ifndef GMXAPI_EXCEPTIONS_H
#define GMXAPI_EXCEPTIONS_H
/*! \defgroup gmxapi_exceptions Exceptions
 *
 * \brief Exceptions thrown by gmxapi components.
 *
 * \ingroup gmxapi
 */
/*! \file
 * \brief Declare exception classes for external API.
 *
 * The core gmxapi library should not throw exceptions to client code. If this happens, please
 * report the bug.
 *
 * The gmxapi header library / template code may use exceptions descended from gmxapi::Exception to
 * simplify client code.
 *
 * In general, errors, warnings, and peculiar circumstances are indicated with the return of empty
 * objects that evaluate false under boolean conversion or with a gmxapi::Status object.
 *
 * \ingroup gmxapi_exceptions
 */

#include <exception>
#include <string>
#include <utility>

namespace gmxapi
{

/*! \brief Base exception for gmxapi library.
 *
 * Exceptions thrown in the gmxapi namespace are descended from gmxapi::Exception
 * or there is a bug.
 *
 * \ingroup gmxapi_exceptions
 */
class Exception : public std::exception
{
    public:
        Exception() = default;
        ~Exception() override = default;
        Exception(const Exception&) = default;
        Exception& operator=(const Exception&) = default;

        Exception(Exception&&) noexcept = default;
        Exception& operator=(Exception&&) noexcept = default;

        const char* what() const noexcept override
        {
            return "Gromacs API error";
        };
};

/*!
 * \brief Basic implementation mix-in for exceptions.
 *
 * Allow exceptions to be defined with minimal syntax when their primary function is
 * to exist as distinct named types.
 *
 * \tparam E the class using this template as a base class.
 *
 * Use in the "curiously recurring template pattern".
 *
 * Example:
 *
 *     class DerivedException : public BasicException<DerivedException>
 *     {
 *          public:
 *              using BasicException::BasicException;
 *     };
 *
 * \note Current implementation only provides constructors and no specialized or dispatched behavior.
 *
 * \ingroup gmxapi_exceptions
 */
template<class E>
class BasicException : public Exception
{
    private:
        std::string what_;
    public:
        BasicException() : BasicException{std::string()}
        {};

        explicit BasicException(std::string&& message) noexcept :
            what_{std::move(message)}
        {};

        explicit BasicException(const char* message)
        {
            what_ = message;
        }

        const char* what() const noexcept override
        {
            return what_.c_str();
        }
};

/*! \brief Behavioral protocol violated.
 *
 * Indicates that a behavioral protocol specified in the API is not being followed. The class
 * throwing this exception expects certain methods to be called in a certain order.
 *
 * If this exception is encountered in client code, the API is being misused or there is a bug.
 * Generally, required behaviors should be implemented in templates or base classes rather than
 * exposing and requiring complete implementation of the protocol in client code.
 *
 * \ingroup gmxapi_exceptions
 */
class ProtocolError : public BasicException<ProtocolError>
{
    public:
        using BasicException<ProtocolError>::BasicException;
};

/*!
 * \brief Intended feature is not implemented.
 *
 * Indicates a bug in the API implementation. Either a version mismatch between the client
 * and library has gone undetected, or the API has purported to offer functionality that does
 * not exist.
 *
 * \ingroup gmxapi_exceptions
 */
class NotImplementedError : public BasicException<NotImplementedError>
{
    public:
        using BasicException<NotImplementedError>::BasicException;
};


/*!
 * \brief Key was not found in a look-up operation
 *
 * The client correctly performed a look-up operation, but the key provided by the client was not
 * found. This is exception is not necessarily an error.
 *
 * \ingroup gmxapi_exeptions
 */
class KeyError : public BasicException<KeyError>
{
    public:
        using BasicException<KeyError>::BasicException;


};

}      // end namespace gmxapi

#endif // header guard
