#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2019, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

"""Abstract base classes for gmxapi Python interfaces.

This module consolidates definitions of some basic interfaces in the gmxapi
Python package. These definitions are evolving and mostly for internal use, but
can be used to check compatibility with the gmxapi implementation details that
are not otherwise fully specified by the API.
"""

import abc
import collections
import typing


class EnsembleDataSource(abc.ABC):
    """A single source of data with ensemble data flow annotations.

    Note that data sources may be Futures.

    Attributes:
        dtype (type): The underlying data type provided by this source.
        source: object or Future of type *dtype*.
        width: ensemble width of this data source handle.

    ..  todo::
        This class should be subsumed into the core gmxapi data model. It is
        currently necessary for some type checking, but will probably disappear
        in future versions.
    """

    def __init__(self, source=None, width=1, dtype=None):
        self.source = source
        self.width = width
        self.dtype = dtype

    @abc.abstractmethod
    def member(self, member: int):
        """Extract a single ensemble member from the ensemble data source."""
        return self.source[member]

    @abc.abstractmethod
    def reset(self):
        """Reset the completion status of this data source.

        Deprecated. This is a workaround until the data subscription model is
        improved. We need to be able to fingerprint data sources robustly, and
        to acquire operation factories from operation handles. In other words,
        a Future will need both to convey its unique recreatable identity as
        well as to be able to rebind to its subscriber(s).

        Used internally to allow graph edges to be reused without rebinding
        operation inputs.
        """
        protocols = ('reset', '_reset')
        for protocol in protocols:
            if hasattr(self.source, protocol):
                getattr(self.source, protocol)()


class NDArray(collections.abc.Sequence):
    """N-dimensional data interface.


    """


# Use SourceTypeVar and ResultTypeVar for static type hints, annotations, and as a parameter to generics.
# Use valid_source_types and valid_result_types for run-time type checking.
ResultTypeVar = typing.TypeVar('ResultTypeVar', *(str, bool, int, float, dict, NDArray))
valid_result_types = ResultTypeVar.__constraints__

SourceTypeVar = typing.TypeVar('SourceTypeVar',
                               *(str, bool, int, float, dict, NDArray, EnsembleDataSource))
valid_source_types = SourceTypeVar.__constraints__
