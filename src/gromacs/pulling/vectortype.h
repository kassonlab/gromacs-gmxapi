//
// Created by Eric Irrgang on 10/13/17.
//

#ifndef GROMACS_VECTORTYPE_H
#define GROMACS_VECTORTYPE_H

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

        vec3() : x{Scalar(0)}, y{Scalar(0)}, z{Scalar(0)} {};
        vec3(const vec3&) = default;
        vec3& operator=(const vec3&) = default;

        /*!
         * \brief Implicit non-narrowing conversion between vec3<>
         * \tparam T
         * \return
         */
        template<typename T>
        operator vec3<T>() { return vec3<T>(x, y, z); }

        /*!
         * \brief Construct a vec3 from arbitrary numeric arguments, but disallow narrowing.
         *
         * \tparam SX numeric type of x
         * \tparam SY numeric type of y
         * \tparam SZ numeric type of z
         * \param X
         * \param Y
         * \param Z
         */
        template<typename SX, typename SY, typename SZ>
        vec3(SX X, SY Y, SZ Z) : x{X}, y{Y}, z{Z} {};
};

//
// Arithmetic
//
// A common idiom in vector math libraries is to overload operator*(), operator/(), and operator%(),
// but in the context of matrices and tensor algebra, it is not unambiguous whether multiplication
// should imply dot product. Let's stick to explicit free functions.
//

template<typename Scalar>
inline vec3<Scalar> operator+(vec3<Scalar> a, vec3<Scalar> b)
{
    return vec3<Scalar>(a.x + b.x, a.y + b.y, a.z + b.z);
};

template<typename Scalar>
inline vec3<Scalar> operator-(vec3<Scalar> a, vec3<Scalar> b)
{
    return vec3<Scalar>(a.x - b.x, a.y - b.y, a.z - b.z);
}

template<typename Scalar>
inline Scalar dot(const vec3<Scalar>& a, const vec3<Scalar>& b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
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

} // end namespace gmx::detail
} // end namespace gmx

#endif //GROMACS_VECTORTYPE_H
