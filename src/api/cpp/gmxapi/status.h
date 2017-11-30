//
// Created by Eric Irrgang on 11/14/17.
//

#ifndef GMXAPI_STATUS_H
#define GMXAPI_STATUS_H
/*! \file
 * \brief Declare the base class for Status returned by API calls.
 *
 * \ingroup gmxapi
 */

#include <memory>

namespace gmxapi{
/*! \brief Container for results of API operations.
 *
 * \internal
 * I'm leaning towards saying this should not be derivable, but that it
 * may contain one or more additional objects, perhaps including exceptions
 * or chains of status / exceptions. Maybe it is a stack. Maybe all
 * API objects should have a Status member that can accumulate Status
 * objects of child objects/operations.
 */
class Status final
{
    public:
        /*!
         * \brief Default constructor.
         */
        Status();
        /*!
         * \brief Copy constructor
         * \param status
         */
        Status(const Status &status);
        /*!
         * \brief Move constructor.
         * \param status
         */
        Status(Status &&status) noexcept;
        /*!
         * \brief Copy assignment operator.
         * \param status
         * \return reference to lhs.
         */
        Status &operator=(const Status &status);
        /*!
         * \brief Move assignment operator.
         * \param status
         * \return reference to lhs
         */
        Status &operator=(Status &&status) noexcept;
        /*!
         * \brief Converting assignment operator.
         * \param success
         * \return reference to lhs
         */
        Status &operator=(bool success);

        /*!
         * \brief Converting constructor.
         * \param success
         */
        explicit Status(bool success);

        /*!
         * \brief non-virtual destructor
         *
         * Do not inherit from this class.
         */
        ~Status();
        /*
         * \brief Check success status.
         *
         * \return true if the operation described was successful.
         */
        bool success() const;
    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
};

} // end namespace gmxapi

#endif //GMXAPI_STATUS_H
