/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */


#include "gmxpre.h"

#include "data.h"

#include "pybind11/pybind11.h"

namespace gmx
{
namespace pyapi
{


LocalTrajDataHandle::LocalTrajDataHandle(shared_ptr<Data3> data) :
    data_(data)
{
}

LocalTrajDataHandle::~LocalTrajDataHandle() {};

/*const shared_ptr<Data3> LocalTrajDataHandle::fetch_data() const
   {
    return data_;
   }*/
shared_ptr<Data3> LocalTrajDataHandle::fetch_data()
{
    return data_;
}

// Export data classes
void export_datatypes(py::module &m)
{
    // Export buffer class that exposes trajectory data arrays
    py::class_< Data3,
                std::shared_ptr<Data3 >
                >(m, "TrajData3", py::buffer_protocol())
    // A buffer interface exported to Python.
    // I'm not sure that def_buffer implies return_value_policy::reference_internal,
    // which implies keep_alive<0,1> to make sure that C++ will keep the
    // object alive as long as the return value is alive, def_buffer does
    // not take return_value_policy arguments. It would be nice if we could
    // instead use the support for Eigen / numpy compatibility.
    // Note that the bindings for std::vector use reference_internal for __getitem__
    // and keep_alive<0,1> for __iter__ but nothing for the def_buffer interface to
    // std::vector. However, it copies data from the py::buffer in the vector_buffer __init__.
    // We should probably use Eigen for data on the C++ side and export directly
    // to numpy arrays instead.
        .def_buffer(
            [](Data3 &data) -> py::buffer_info
            {
                return py::buffer_info(
                        data.data(),                                /* Pointer to buffer */
                        sizeof(real),                               /* Size of one scalar */
                        py::format_descriptor<real>::format(),      /* Python struct-style format descriptor */
                        2,                                          /* Number of dimensions */
                        { data.N(), data.dim() },                   /* Python buffer dimensions */
                        { sizeof(real) * data.dim(), sizeof(real) } /* Strides (in bytes) for each index in C++ */
                        );
            }
            )
    // Accept buffers from Python.
    // We should whether and how to safely perform a no-copy construction from a
    // buffer when TrajDataArray can have multiple references on the C++ side.
    // If the Python buffer views are all closed and there are no more
    // Python references to the object, then any remaining C++ references
    // to the object will have their data become invalid if the buffer object
    // is allowed to be released. I need to look into this more, but if pybind11 doesn't
    // already do it, we can keep the source of the buffer alive by using
    // keep_alive<1,3> in the init method to keep the py::buffer argument
    // alive as long as the TrajDataArray. return_value_policy::reference
    // prevents Python from taking ownership (C++ manages the lifetime).
    // That may interfere with our ability to release the global interpreter
    // lock. Instead, we should probably only attempt non-temporary objects
    // with non-copy construction on memory that is already owned by the API and
    // not the Python interpreter or retain handles to Python objects (providing
    // a buffer interface) that are only exposed in controlled and limited scope.
    // If we want to set
    // data in TrajDataArray objects with minimal copies, we allocate our own
    // memory and implement
    // element access methods to write directly from Python to the managed
    // array.
    // TO DO: only accept dense numpy arrays with array_t arguments.
    // py::array_t<real, py::array::c_style | py::array::forcecast> data
        .def("__init__",
             [](Data3 &data, py::buffer b)
             {
                 /* Request a buffer descriptor from Python */
                 py::buffer_info info = b.request();

                 /* Some sanity checks ... */
                 if (info.format != py::format_descriptor<real>::format())
                 {
                     throw std::runtime_error("Incompatible format: expected a array of type real!");
                 }
                 ;
                 if (info.ndim != 2 || info.shape[0] != 3)
                 {
                     throw std::runtime_error("Incompatible buffer dimension!");
                 }
                 ;

                 // Construct in place
                 // It is important that the reference count of the buffer object b
                 // should be incremented to prevent Python garbage collection from
                 // deallocating the memory in the TrajDataArray object. I assume
                 // pybind11 takes care of that.
                 new (&data)Data3(static_cast<real *>(info.ptr), info.shape[0]);
             },
             py::keep_alive<1, 3>() // keep py::buffer b alive while *this is alive
             )
    // Inspect...
        .def("__repr__", [](const Data3 &t)
             {
                 std::stringstream repr;
                 //std::string repr{"tah dah!"};//
                 repr << t.N() << "x" << t.dim() << " array of trajectory data of type 'real'\n";
                 // Ugh...
                 for (size_t i = 0; i < t.N(); ++i)
                 {
                     repr << t[i][0] << "\t" << t[i][1] << "\t" << t[i][2] << "\n";
                 }
                 return repr.str();
             }
             )
        .def_property_readonly("N", &Data3::N, "number of elements")
    /* Needs better TrajDataArray
       .def("__iter__", [](const Data3 &s)
         {
             return py::make_iterator(s.begin(), s.end());
         },
         py::keep_alive<0, 1>() // Essential: keep object alive while iterator exists
        )
     */
    // Use generator or make a list instead...
    //.def("__getitem__", [](const Data3 &s, py::slice slice) -> Data3* {})
    ;

    // Export an API object to serve as a light-weight handle.
    py::class_< Data3Handle > (m, "GmxData3Base");
    py::class_< LocalTrajDataHandle, Data3Handle > (m, "GmxData3")
        .def("extract", &LocalTrajDataHandle::fetch_data, "Extract API object to Python interpreter");
    ;

}

} // end namespace pyapi
} // end namespace gmx
