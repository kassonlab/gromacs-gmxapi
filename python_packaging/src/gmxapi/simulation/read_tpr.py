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

"""read_tpr operation module

Provides implementation classes and user interface for gmxapi.read_tpr.
"""

import inspect
import typing

import gmxapi
import gmxapi.fileio
import gmxapi.operation as _op


class OutputDataProxy(_op.DataProxyBase, descriptors={'parameters': _op.OutputDescriptor('parameters', dict)}):
    """Implement the 'output' attribute of ReadTpr operations."""


class Publisher(gmxapi.operation.Publisher):
    """Implement the publishing data descriptor for ReadTpr parameters output."""
    def __set__(self, instance: 'PublishingDataProxy', value):
        super().__set__(instance, value)


class PublishingDataProxy(_op.DataProxyBase,
                          descriptors={
                              'parameters': Publisher('parameters', dict)}
                          ):
    """Manage output resource updates for ReadTpr operation."""


# Note: we borrow the implementation from operation.ResourceManager for now,
# but in the future we want the implementations to either be decoupled or
# for implementations in a given context to be coupled to details that are clearly
# and explicitly related to that context. Right now, operation.ResourceManager
# is tied to the implementation of Contexts in gmxapi.operation, but that is not
# sufficiently clear and explicit.
class ResourceManager(gmxapi.operation.ResourceManager):
    """Manage resources for the ReadTpr operation in the gmxapi.operation contexts.

    Extends gmxapi.operation.ResourceManager to tolerate non-standard data payloads.
    Futures managed by this resource manager may contain additional attributes.
    """
    def future(self, name: str, description: _op.ResultDescription):
        tpr_future = super().future(name=name, description=description)
        return tpr_future

    def data(self) -> OutputDataProxy:
        return OutputDataProxy(self)


class Resources(object):
    """Input and output run-time resources for a ReadTpr operation."""
    def __init__(self, tpr_filename, publisher: PublishingDataProxy):
        self.tpr_object = gmxapi.fileio.TprFile(filename=tpr_filename, mode='r')
        self.output = publisher

    def filename(self):
        return self.tpr_object.filename

    @staticmethod
    def create(input, output) -> 'Resources':
        """Factory function to create new resource instances."""
        filename = input.kwargs['filename']
        return Resources(tpr_filename=filename, publisher=output)


class ReadTprDetails(_op.OperationDetailsBase):
    """Produce simulation input data sources from a TPR file.

    Outputs:
        parameters : simulation parameters
        topology : molecular force field data
        coordinates : phase space coordinates
        simulation_state : integrator state and other checkpoint data

    The output is suitable as input for operations such as modify_input or mdrun.
    """
    _output = _op.OutputCollectionDescription(parameters=dict)
    _input = _op.InputCollectionDescription(
        [('filename', inspect.Parameter('filename',
                                        inspect.Parameter.POSITIONAL_OR_KEYWORD,
                                        annotation=str))])

    def __call__(self, resources: Resources):
        # Operation implementation in the gmxapi.operation module context.
        with resources.tpr_object as fh:
            params = fh._tprFileHandle.params().extract()
            resources.output.parameters = params

    @classmethod
    def signature(cls) -> _op.InputCollectionDescription:
        return cls._input

    @classmethod
    def resource_director(cls, *, input, output) -> Resources:
        resources = Resources.create(input, output)
        return resources

    def output_description(self) -> _op.OutputCollectionDescription:
        return self._output

    def publishing_data_proxy(self, *, instance: ResourceManager, client_id) -> PublishingDataProxy:
        return PublishingDataProxy(instance=instance, client_id=client_id)

    def output_data_proxy(self, instance: ResourceManager) -> OutputDataProxy:
        return OutputDataProxy(instance=instance)

    @classmethod
    def make_uid(cls, input) -> str:
        # TODO: input.fingerprint()
        # TODO: handle edge width
        filename = input.source_collection['filename']
        return 'read_tpr_' + filename

    @classmethod
    def operation_director(cls, *args, context: _op.Context, label=None, **kwargs) -> _op.AbstractOperation:
        return super().operation_director(*args, context=context, label=label, **kwargs)


class ReadTpr(_op.AbstractOperation):
    def __init__(self, resource_manager: ResourceManager):
        self.__resource_manager = resource_manager

    def run(self):
        self.__resource_manager.update_output()

    @property
    def output(self) -> OutputDataProxy:
        return self.__resource_manager.data()


def read_tpr(filename=None):
    # TODO: Use NodeBuilder protocol to add the operation instance to the Context.
    # context = current_context()
    # return ReadTpr.operation_director(filename=filename, context=context)
    operation = ReadTprDetails()
    sources = _op.DataSourceCollection()
    sources['filename'] = filename
    input_sink = _op.SinkTerminal(ReadTprDetails.signature())
    input_sink.update(sources)
    edge = _op.DataEdge(sources, input_sink)

    resource_manager = ResourceManager(operation=operation, source=edge)
    return ReadTpr(resource_manager)

# # TODO: fix make_operation
# # read_tpr = make_operation(fileio._SimulationInput, output={'parameters': dict})
# @function_wrapper(output={'parameters': dict, 'topology': str, 'coordinates': str, 'simulation_state': str})
# def read_tpr(filename: str = '', output=None):
#     """Get simulation input sources from a TPR file.
#
#     Outputs:
#         parameters : MDP simulation parameters
#         coordinates : atom (or CG particle) coordinates (not yet implemented)
#         simulation_state : simulation internal state (checkpoint data) (not yet implemented)
#         topology : molecular force field data (not yet implemented)
#     """
#     import gmxapi.fileio
#     sim_input = gmxapi.fileio.read_tpr(filename)
#     output._internal = sim_input
#     output.parameters = sim_input.parameters.extract()
#     output.topology = filename
#     output.coordinates = filename
#     output.simulation_state = filename
