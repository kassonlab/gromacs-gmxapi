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

Refer to [PEP 484](https://www.python.org/dev/peps/pep-0484) and to
[PEP 526](https://www.python.org/dev/peps/pep-0526) for Python background.

Developer Note:
    Type checking fails when accessing attributes of generic classes (see PEP 484).
    This affects some scoping choices.

Developer Note:
    It is worth noting some details regarding type erasure with Python generics.
    Generic type parameters (class subscripts) allow for annotation of dynamic
    types that are bound when the generic is instantiated. The subscript is not
    part of the run time class. This means that we can use subscripted generic
    classes for static type checking, but use for any other purpose is discouraged.
    Thus, type variable parameters to generics have less meaning that do template
    parameters in C++, and are orthogonal to subclassing. Note, though, that a
    subclass is not generic if the generic parameters are bound (explicitly
    with type subscripts or implicitly by omission (implying `typing.Any`).

    Also, type parameters
    (or types of parameters) to functions can be used to dispatch a generic
    function to a specific function.

    In other words: keep in mind the
    orthogonality of generic classes and base classes, and recognize that
    composed objects mimicking C++ template specializations are not distinguishable
    at the class level.

"""
# Note that the Python typing module defines generic classes in terms of abstract
# base classes defined in other modules (namely `collections`), but without
# actual inheritance. The ABCs are not intended to be instantiated, and the
# generics are very mangled objects that cannot be instantiated. However,
# user-defined subclasses of either may be instantiated.
#
# This restriction is probably due somewhat to implementation constraints in the
# typing module, but it represents a Separation of Concerns that we should
# consider borrowing in our model. Practically, this can mean using an abstract
# for run time checking and a generic for static type checking.
#
# There may not be a compelling reason to
# rigorously separate ABC and generic type, or to disallow instantiating a
# generic unless it is also an abstract. Note the distinction, though, between
# abstract generics and fully implemented generics. NDArray
# and Future are likely to be examples of fully implemented generics, while
# Context and various interface types are likely to have abstract generics.
#
# Use abstract base classes to define interfaces and relationships between
# interfaces. Use `typing` module machinery for static type checking. Use the
# presence or absence of an expected interface, or exception handling, for run
# time type checking. Use `isinstance` and `issubclass` checking against
# non-generic abstract base classes when the valid interface is complicated or
# unknown in the caller's context.
import typing
from abc import ABC, ABCMeta, abstractmethod
import collections
from typing import TypeVar, Generic, NewType, Type, Callable


class EnsembleDataSource(ABC):
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

    @abstractmethod
    def member(self, member: int):
        """Extract a single ensemble member from the ensemble data source."""
        return self.source[member]

    @abstractmethod
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


class NDArray(collections.abc.Sequence, ABC):
    """N-dimensional data interface.


    """


# Use SourceTypeVar and ResultTypeVar for static type hints, annotations, and as a parameter to generics.
# Use valid_source_types and valid_result_types for run-time type checking.
ResultTypeVar = TypeVar('ResultTypeVar', *(str, bool, int, float, dict, NDArray))
valid_result_types = ResultTypeVar.__constraints__

SourceTypeVar = TypeVar('SourceTypeVar',
                        *(str, bool, int, float, dict, NDArray, EnsembleDataSource))
valid_source_types = SourceTypeVar.__constraints__

# Place holder for type annotations of Context objects.
# TODO: Expand to support some static type checking.
_Context = TypeVar('_Context')
# Type variable that binds to subclasses of the forward referenced OperationImplentation ABC.
_Op = TypeVar('_Op', bound='OperationImplementation')


class Future(ABC):
    """Generic result data."""
    @property
    @abstractmethod
    def dtype(self) -> Type[typing.Union[valid_result_types]]:
        ...

    @abstractmethod
    def result(self) -> typing.Union[valid_result_types]:
        ...


class OutputDataProxy(ABC):
    """A collection of Operation outputs.

    This abstract base class describes the interface to the output of an
    operation / work node.
    """
    # TODO: Specification.
    # Currently, the common aspect of OutputDataProxy is that a class has public
    # attributes that are exclusively OutputDescriptor instances, meaning that
    # getattr(instance, attr) returns a Future object. However, there are several
    # ways to implement getattr, and this sort of check in an ABC does not appear
    # to be common.
    #
    # We might choose to assert that all public attributes must be compatible data
    # descriptors, in conjunction with defining a more specific OutputDataProxy
    # metaclass, but this makes for a dubiously growing chain of data descriptors
    # we use for the output access.
    #
    # The data model might be cleaner if we move to something more
    # like a Collection or Mapping with more conventional getters, but we would
    # lose the ability to check type on individual elements. (However, the
    # typing of return values is not normally the defining aspect of an ABC.)
    #
    # Another alternative is to collapse the contents of the `output` attribute
    # into the Operation handle type, strongly define all handle types (so that
    # the type checker can identify the presence of attributes), and rely only
    # on type checking at the level of the data descriptors. (Dynamically defined
    # OutputDataProxy classes are the execption, rather than the rule.)
    #
    # We will need to consider the details of type checkers and syntax inspection
    # tools, like Jedi and MyPy, to make design choices that maximize API usability
    # and discoverability.


class AbstractOperationReference(ABC):
    """Client interface to an element of computational work already configured.

    An "instance" of an operation is assumed to be a node in a computational
    work graph, owned and managed by a Context. This class describes the
    interface of the reference held by a client once the node exists.
    """

    @abstractmethod
    def run(self):
        """Assert execution of an operation.

        After calling run(), the operation results are guaranteed to be available
        in the local context.
        """

    @property
    @abstractmethod
    def output(self) -> OutputDataProxy:
        """Get a proxy collection to the output of the operation.

        Developer note: The 'output' property exists to isolate the namespace of
        output data from other operation handle attributes and we should consider
        whether it is actually necessary or helpful. To facilitate its possible
        future removal, do not enrich its interface beyond that of a collection
        of OutputDescriptor attributes.
        """
        ...


class OperationReference(AbstractOperationReference, Generic[_Op]):
    """Object with an OperationReference interface.

    Generic version of AbstractOperationReference, parameterized by operation
    implementation.
    """


class NodeBuilder(ABC):
    """Add an element of computational work to be managed by a gmxapi Context."""

    @abstractmethod
    def build(self) -> AbstractOperationReference:
        """Finalize the creation of the operation instance and get a reference."""
        ...

    @abstractmethod
    def add_input(self, name: str, source):
        """Attach a client-provided data source to the named input."""
        ...

    @abstractmethod
    def add_resource_factory(self, factory: Callable):
        """Register a resource factory for the operation run-time resources.

        The factory will be called within the Context
        """
        # The factory function takes input in the form the Context will provide it
        # and produces a resource object that will be passed to the callable that
        # implements the operation.
        assert callable(factory)
        ...


class Context(ABC):
    """API Context.

    All gmxapi data and operations are owned by a Context instance. The Context
    manages the details of how work is run and how data is managed.

    This abstract base class (ABC) defines the required interface of a Context
    implementation. Client code should use this ABC for type hints. Concrete
    implementations may, *but are not required*, to subclass from this ABC to
    help enforce compatibility.
    """
    @abstractmethod
    def node_builder(self, *, operation, label=None) -> NodeBuilder:
        """Get a builder for a new work graph node.

        Nodes are elements of computational work, with resources and execution
        managed by the Context. The Context handles parallelism resources, data
        placement, work scheduling, and data flow / execution dependencies.

        This method is used by Operation director code and helper functions to
        add work to the graph.
        """
        ...


class OperationDirector(Generic[_Op, _Context]):
    """Generic abstract operation director.

    An operation director is instantiated for a specific operation and context
    by a dispatching factory to add a computational element to the work managed
    by the context.

    Note:
        This is a generic class, as defined with the Python `typing` module.
        If used as a base class, re-expression of the TypeVar parameters in class
        subscripts will cause the derived class to be generic unless regular
        classes (non-TypeVar) are given. Omission of the subscripts causes `Any`
        to be bound, which also results in a non-generic subclass.
    """

    # TODO: Annotate `input`, whose validity is determined by both context and operation.
    @abstractmethod
    def __call__(self, input):
        """Add an element of work (node) and return a handle to the client."""
        ...

    def handle_type(self) -> Type[OperationReference[_Op]]:
        """Get the class used for operation references in this Context."""
        ...


class OperationImplementation(ABC):
    """Essential interface of an Operation implementation.

    Describe the essential features of an Operation that can be registered with
    gmxapi to support building and executing work graphs in gmxapi compatible
    execution contexts.
    """
    # The executable part of an operation consumes a distinct resource type.
    # The resource type may be opaque, because it is created through a factory
    # and only used in passing to a function object.
    ResourceType = NewType('ResourceType', object)

    # TODO: Consider a data descriptor and metaclass to validate the name and namespace.
    @classmethod
    @abstractmethod
    def name(cls) -> str:
        """The name of the operation.

        Generally, this corresponds to a callable attribute of a Python module
        (named by namespace()) that acts as a factory for new operation instances.
        It is also used by Context implementations to locate code supporting
        the operation.
        """
        ...

    # TODO: Consider a data descriptor and metaclass to validate the name and namespace.
    @classmethod
    @abstractmethod
    def namespace(cls) -> str:
        """The namespace of the operation.

        Generally, the namespace corresponds to a Python module importable in
        the execution environment.
        """

    # Note that this indicates that the ABC requires subclasses to provide a generic function,
    # which _is_ the factory behavior we are trying to specify.
    # TODO: Confirm that this type-checks correctly.
    # We may need to look more at how the type checking is implemented to see how
    # to do this right.
    @classmethod
    @abstractmethod
    def director(cls: Type[_Op], context: _Context) -> OperationDirector[_Op, _Context]:
        """Factory to get an OperationDirector appropriate for the context."""
        ...
