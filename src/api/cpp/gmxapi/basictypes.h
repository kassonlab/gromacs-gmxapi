#ifndef GMXAPI_BASIC_TYPES_H
#define GMXAPI_BASIC_TYPES_H
/// \cond DEV
/*! \file
 * \brief Basic data structures, formats, and access types.
 *
 * Data structures could be templated on fundamental type, memory layout, and access
 * type. Sparse arrays or other non-contiguous arrays could be possible,
 * and would not be convertible to raw pointers without copy, while
 * simple arrays are convertible but have associated layout information.
 *
 * Important layout information includes:
 * * fast indexing direction (row major or column major)
 * * row stride
 * * column stride
 * * nominal width
 * * offset of first element
 *
 * It may be necessary to treat two distinct use cases of light-weight
 * handle acquisition for performance-sensitive code
 */

#include <functional>
#include <memory>
#include <utility>

namespace gmxapi
{

// enum class Access : int
// {
//     READ_ONLY,       // Read-only handle points to data state at time of call.
//     // Class owning data may update contents at a new address during the lifetime
//     // of the handle by implementing copy-on-write, but they may instead default to
//     // locked read semantics.
//     LOCKED_WRITE,    // Zero or one write handles can exist at a time.
//     COPY_ON_WRITE,   // Like a READ_ONLY handle, but allows non-const
//                      // type by guaranteeing copy-on-write.
//     LOCKED_READ      // Handle is guaranteed to be current for lifetime.
//                      // Prevents write handles.
//     // Additional access modes may be necessary. Write-only or sequential-only
//     // handles may be necessary for streams, but it may be better to have
//     // separate a separate streams class.
// };

/*!
 * Acquiring a handle does not guarantee that data is immediately accessible
 * in contiguous memory, and writes to a handle are not guaranteed to be visible
 * outside of the current scope until the handle is released. sync() and flush()
 * calls may be appropriate.
 *
 *
 */
// Is there a sensible use for access modes as types to make templating easier?
// Alternatively, skip the template parameter idea and declare a class for
// each access mode. (Probably loses some compiler optimizations.)
template<class T, class Access> class ManagedHandle
{
    public:
        typedef Access mode;
        typename T::value_type &operator[](typename T::size_type i);
        const typename T::value_type &operator[](typename T::size_type i) const;

//    static const Access mode{a};

        // Provide a means to find the managing object from which to request a new handle, as well as to keep the manager alive.
        //const std::shared_ptr<Managed<T>> ref_;
    private:
        T data_;
};

template<class T> class read_only
{
    public:
        using callback_t = std::function<void()>;

        read_only()                              = delete;
        read_only(read_only &&handle)            = default;
        read_only &operator=(read_only &&handle) = default;

        // Join shared ownership.
        explicit read_only(const std::shared_ptr<T> ptr) : data_ {ptr}
        {};
        // Copy data
        explicit read_only(const T &data) : data_ {std::make_shared<T>(data)}
        {};
        // Take ownership of data
        explicit read_only(T &&data) : data_ {std::shared_ptr<T>(std::move(data))}
        {};

        //
        read_only(std::shared_ptr<T> ptr, std::shared_ptr<callback_t> releaser) :
            data_ {ptr},
        releaser_ {releaser}
        {};

        ~read_only()
        {
            // If the functor that was bound at creation still exists, make the
            // provided call to allow clever book-keeping.
            if (auto release = releaser_.lock())
            {
                (*release)();
            }
        };

        typename T::value_type &operator[](typename T::size_type i);
        const typename T::value_type &operator[](typename T::size_type i) const;
    private:
        const std::shared_ptr<T>  data_;
        std::weak_ptr<callback_t> releaser_;
        //std::weak_ptr<Managed<T>> ref_;
};

template<class T> class locked_write : public ManagedHandle < T, locked_write < T>>
{
};
template<class T> class copy_on_write : public ManagedHandle < T, copy_on_write < T>>
{
};
template<class T> class locked_read : public ManagedHandle < T, locked_read < T>>
{
};

// Partially specialize for access modes.
template<class T>
typename T::value_type &read_only<T>::operator[](typename T::size_type i)
{
// Assume T with value_type and size_type has a subscript operator...
    return (*data_)[i];
}

/*! \brief Provide different sharing behaviors for data.
 *
 * Classes providing managed data handles inherit this when consumers need multiple
 * sharing semantics for extending the lifetime of data, allowing
 * optimized storage, and interacting with data state.
 *
 * Managed objects can be "move"d to transfer ownership, but copying
 * must be explicit to avoid accidental misuse. Managed objects
 * are thread-safe process-local data.
 *
 * It is difficult to predict how data will be used, how it will be passed
 * to other code, and what are the relative lifetimes of consuming objects.
 * Some consumers might need to extend data lifetimes. This should help,
 * as well as providing sanity checking to avoid ambiguities of data
 * state.
 *
 * Managed and ManagedHandle use the Template Method behavioral pattern
 * and/or the Curiously Recurring Template Pattern to achieve performance and
 * reduce dependencies after compile time, as well as to enforce certain behavior in derived classes.
 *
 * Objects should be constructed with initialization to avoid resource allocation
 * surprises or wasteful reallocations later, but default construction is still
 * possible to allow initialization outside of the owning scope.
 *
 * Derived classes do not have a proper "is a" substitional property, but rather
 * are "implemented in terms of".
 *
 * TBD:
 * Does a handle need to extend the life of the manager, or just the data? I.e.
 * should a handle hold a weak_ptr or shared_ptr to the manager?
 * I think the Handle and Manager do not need to be so tightly coupled. particularly
 * if a Handle can create a new Manager from a write handle (or const Manager
 * from a read handle). For this to work, data needs to be managed with a
 * reference-counted smart pointer.
 *
 */
// TODO: check syntax to make sure we do not allow T to be a reference.
template<class T>
class Managed
{
    public:
        // The data can be initialized at construction, but afterwards requires a
        // write handle to modify.
        template<class ... Targs> Managed(Targs && ... args) :
            data_ {std::forward<Targs>(args) ...}
        {};

        Managed(std::initializer_list<typename T::value_type> elements) :
            data_ {std::make_shared<T>(elements)}
        {};

        // Manager can take over ownership of moveable objects or initialize from
        // copies.
        Managed(T &&data) : data_ {std::move(data)}
        {};

        // Manager can be moved to transfer ownership.
        Managed(Managed<T> &&) = default;

        // Manager should not be implicitly copied to avoid misuse.
        Managed(const Managed<T> &) = delete;

        ~Managed() = default;

        read_only<T> handle();
        //template<class Access> Access handle();
        //TODO: handle const managers
        //template<Access a> ManagedHandle<T, a> handle() const;

        // template<>
        // ManagedHandle<T, std::enable_if</*write*/>> handle() const
        // {
        //     // Can't request a write handle from a const manager.
        //     throw APIError;
        // }

        // Temporary raw accessor for initial testing.
        T* get() {return data_.get(); };
        const T* get() const {return data_.get(); };

    private:
        std::shared_ptr<T> data_;
        // write handles and locked reads on const pointers must prevent new write handles.
        //mutable std::mutex writeLock_;
        // Allow multiple locked read handles.
        //mutable std::mutex readLock_;
        // We must track the number of handles to allow optimal memory management.
        //mutable unsigned int nHandles_;
};

template<class T>
read_only<T> Managed<T>::handle()
{
    return read_only<T>(*data_);
}

// // single vectors are simple compact structs that should be used for
// // convenience where they will be optimized away by the compiler.
// template<typename Scalar>
// struct Vector3D
// {
//     Scalar x;
//     Scalar y;
//     Scalar z;
//
//     // Default constructor guarantees zero-initialization
//     Vector3D() :
//         x{Scalar(0)},
//         y{Scalar(0)},
//         z{Scalar(0)}
//     {};
//
//     // Element-wise construction avoids accidental convesion or narrowing.
//     explicit Vector3D(const Scalar& a, const Scalar& b, const Scalar& c) :
//         x{a},
//         y{b},
//         z{c}
//     {};
//
//     // Initializer list constructor
//
//     ~Vector3D() = default;
//     explicit Vector3D(const Vector3D<Scalar>&) = default;
// };
//
// /*! \brief NxD array data.
//  *
//  * Consider Eigen to avoid reinventing the wheel.
//  */
//
//
// /// A compact
// template<typename Scalar>
// std::array<Scalar, 3>::operator Vector3D<Scalar, compact>()
// {
//     return std::array<Scalar, 3>{vec3::x, vec3::y, vec3::z};
// }

}      // end namespace gmxapi
/// \endcond
#endif // header guard
