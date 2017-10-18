//
// Created by Eric Irrgang on 10/13/17.
//

#ifndef GROMACS_VECTORTYPE_H
#define GROMACS_VECTORTYPE_H

/*! \file
 * \brief Template header for 3D vector types and operations.
 *
 * Two reasons:
 *
 * 1. Make data types and precision explicit and unambiguous.
 * 2. Provide an abstraction from storage method.
 *
 * These types should map easily to float3 (or float4) as in CUDA and other libraries,
 * as well as to arrays or even non-contiguous structures, at least insofar as
 * the compiler should be able to optimize away copies.
 *
 * Along these lines, the structures are intended to be short-lived handles
 * for convenience and strong typing of operations. Arrays of vec3 should not
 * be necessary and are probably not desirable, at least across C++ translation units.
 *
 */

// TODO: It's worth taking a look at how compilers handle iterative extraction of vec3 from Nx3 data in practice.

#include <cmath>
#include "gromacs/math/vectypes.h"

namespace gmx
{
namespace detail
{

/*!
 * \brief 3-dimensional vector types.
 *
 * Provide strongly-typed vector type for unambiguous operations. Maps easily to CUDA float3.
 *
 * Conversions to other types with different semantics should be explicit and avoid accidental loss of precision.
 * \tparam Scalar
 */
template<typename Scalar>
class vec3
{
    public:
        Scalar x;
        Scalar y;
        Scalar z;

        /*!
         * \brief Require type matching for direct construction.
         * \param X
         * \param Y
         * \param Z
         */
        constexpr explicit vec3(const Scalar& X, const Scalar& Y, const Scalar& Z) : x{X}, y{Y}, z{Z} {};

        vec3() : vec3{Scalar(0), Scalar(0), Scalar(0)} {};
        vec3(const vec3&) = default;
        vec3& operator=(const vec3&) = default;
        vec3& operator=(vec3&&) noexcept = default;

        /*!
         * \brief Implicit non-narrowing conversion between vec3<>
         * \tparam T
         * \return
         */
        template<typename T>
        explicit operator vec3<T>() { return vec3<T>(x, y, z); }

        template<typename T>
        explicit vec3(const T& a) : vec3(a.x, a.y, a.z) {};
};

//
// Arithmetic
//
// A common idiom in vector math libraries is to overload operator*(), operator/(), and operator%(),
// but in the context of matrices and tensor algebra, it is not unambiguous whether multiplication
// should imply dot product. Let's stick to explicit free functions.
//

/*!
 * \brief Unary negation operator
 * \tparam Scalar underlying vector type
 * \param a input vector
 * \return (-a.x, -a.y, -a.z)
 */
template<typename Scalar>
inline vec3<Scalar> operator-(const vec3<Scalar>& a)
{
    return vec3<Scalar>(-a.x, -a.y, -a.z);
}

template<typename Scalar>
inline vec3<Scalar> operator+(const vec3<Scalar>& a, const vec3<Scalar>& b)
{
    return vec3<Scalar>(a.x + b.x, a.y + b.y, a.z + b.z);
};

template<typename Scalar>
inline vec3<Scalar> operator-(const vec3<Scalar>& a, const vec3<Scalar>& b)
{
    return vec3<Scalar>(a.x - b.x, a.y - b.y, a.z - b.z);
}


/*!
 * \brief Multiply vector by scalar
 *
 * \tparam Scalar vector type
 * \param a input vector
 * \param s input scalar
 * \return (a.x, a.y, a.z) * s
 *
 * Note that input scalar may be implicitly narrowed if higher precision than input vector.
 */
template<typename T1, typename T2>
inline vec3<T1> operator*(const vec3<T1>& a, const T2& s)
{
    return vec3<T1>(T1(a.x * s), T1(a.y * s), T1(a.z * s));
}
template<typename T1, typename T2>
inline vec3<T1> operator*(const T2& s, const vec3<T1>& a)
{
    return a * s;
};

/*!
 * \brief Vector division by scalar
 * \tparam Scalar underlying vector type
 * \param a input vector
 * \param b input scalar
 * \return
 *
 * Note that input scalar may be implicitly narrowed if higher precision than input vector.
 */
template<typename T1, typename T2>
inline vec3<T1> operator/(const vec3<T1>& a, const T2& s)
{
    assert(s != T2(0.0));
    const T2 inv_s{T2(1.0)/s};
    return vec3<T1>(T1(a.x * inv_s), T1(a.y * inv_s), T1(a.z * inv_s));
}

/*!
 * \brief Scalar vector product
 * \tparam Scalar underlying vector type
 * \param a first vector
 * \param b second vector
 * \return (a.x * b.x) + (a.y * b.y) + (a.z * b.z)
 */
template<typename Scalar>
inline Scalar dot(const vec3<Scalar>& a, const vec3<Scalar>& b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

/*!
 * \brief Norm or magnitude of vector
 */
template<typename Scalar, typename Tout = double>
inline Tout norm(const vec3<Scalar>& a)
{
    return sqrt(dot(vec3<Tout>(a), vec3<Tout>(a)));
}

//
// Comparisons
//

/*!
 * \brief Equality comparison operator
 * \tparam Scalar
 * \param a
 * \param b
 * \return true if all elements of vectors are arithmetically equal.
 */
template<typename S1, typename S2>
bool operator==(const vec3<S1>& a, const vec3<S2>& b)
{
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

//
// Conversions
//

template<typename Scalar>
inline gmx::RVec as_Rvec(const vec3<Scalar>& v)
{
    return {v.x, v.y, v.z};
}

//
// Helpers
//

/*!
 * \brief Flexibly produce vector of a given type.
 *
 * \tparam Scalar underlying output vector type
 * \tparam T1 x input type
 * \tparam T2 y input type
 * \tparam T3 z input type
 * \param x
 * \param y
 * \param z
 * \return vector (x, y, z) of type vec3<Scalar>
 *
 * Helper function allows narrowing and mismatched types. Constructs a vec3<Scalar> from any x, y, and z that can be
 * implicitly converted to type Scalar.
 */
template<typename Scalar, typename T1, typename T2, typename T3>
inline constexpr vec3<Scalar> make_vec3(T1 x, T2 y, T3 z)
{
    return vec3<Scalar>(Scalar(x), Scalar(y), Scalar(z));
};


} // end namespace gmx::detail
} // end namespace gmx

#endif //GROMACS_VECTORTYPE_H
