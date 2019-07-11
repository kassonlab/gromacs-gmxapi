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

"""GROMACS simulation subpackage for gmxapi.

Contains implementation details for launching and manipulating simulation work.

The mdrun operation (in its first draft) conforms to the user-level API, but does
not use the Python Context resource manager. It uses either the legacy 0.0.7
Context or its own Context, also implemented in this module.
"""

__all__ = ['mdrun', 'read_tpr']

import inspect
import typing

from gmxapi import exceptions
from gmxapi.operation import SourceResource, DataProxyBase, OutputDescriptor, ResultDescription, OutputData, \
    AbstractOperation, OperationDetailsBase
# The following imports are not marked as public API.
# TODO: Resolve public API and exposure.
from gmxapi.operation import current_context, function_wrapper, InputCollectionDescription, OutputCollectionDescription
from gmxapi.simulation import workflow


def mdrun_dispatcher(context, *, input, label: str = None, **kwargs) -> 'MDRun':
    """Dispatch to an appropriate director based on the context and input.

    Runs appropriate director code to set up an operation, returning a handle to
    a simulation operation node.
    """
    from gmxapi.simulation.context import get_context
    if label is not None:
        raise exceptions.NotImplementedError('sorry... no labels yet')
    try:
        legacy_context = get_context(work=input)
    except Exception:
        legacy_context = None

    def run_session():
        with legacy_context as session:
            session.run()
        return True

    if context is not None and context is legacy_context:
        helper = function_wrapper(output={'data': bool})(run_session)
        return helper(**kwargs)
    else:
        raise exceptions.ValueError('Could not dispatch MD input {} with context {}'.format(input, legacy_context))


# class MDRunDirector(object):
#     """Dispatch instantiation of implementation objects, informed by context."""


class MDRunResourceDirector(object):
    """Provide factory for MDRun operation run-time resources.

    Collaborates with resource manager to separate details of resource management
    or execution context from the operation implementation.
    """


class MDRunImplementation(OperationDetailsBase):
    """Provide the mdrun operation in the gmxapi python package context.
    """


class MDRunResourceManager(SourceResource):
    """Resource management gateway for mdrun outputs to Python data flow."""

    class MDRunOutputData(DataProxyBase):
        """Provide the "output" attribute for mdrun operations.

        'coordinates' is a mapping for microstate phase space coordinates. Contains
        'positions' for particle positions (Nx3).
        """
        # First draft: no output.
        coordinates = OutputDescriptor('coordinates', dtype=dict)
        parameters = OutputDescriptor('parameters', dtype=dict)

    coordinates_description = ResultDescription(dtype=dict)
    parameters_description = ResultDescription(dtype=dict)

    def __init__(self):
        self._data = OutputData('coordinates', description=self.coordinates_description)

    def data(self) -> MDRunOutputData:
        return self.MDRunOutputData(instance=self)

    def is_done(self, name: str) -> bool:
        return self._data.done

    def get(self, name: str) -> OutputData:
        if name != self._data.name:
            raise exceptions.ValueError('Unknown resource requested.')
        return self._data

    def update_output(self):
        # We need to make sure that the simulation runs no more than once within
        # the resources available to the client. To that end, we
        pass

    def reset(self):
        pass

    def width(self) -> int:
        return self._data._description.width


class MDRun(AbstractOperation):
    """Handle to an MD simulation.

    In the default Context, can take filename argument for TPR file
    and produce filename output for trajectory output.

    In the context of a compatible operation, can provide additional non-standard
    outputs. In the context of a ReadTpr or ModifyInput operation, consumes
    structured input that can be raw TPR file contents with possibly modified parameters.

    Provides additional non-standard methods for the Operation.

    Attach an MD plugin...

    .. uml::


    This class is not part of the public interface.
    Client gets instances of this through gmx.mdrun() factory.
    """

    def __init__(self):
        """Create a gmxapi graph node with extra GROMACS library bindings.

        C++ bindings are not standardized yet, and the input to initialize an
        MD operation is not very stable.
        """

    def run(self):
        # The instruction to run is proxied to the resource manager.
        pass

    @property
    def output(self) -> DataProxyBase:
        """Provide MD outputs.

        E.g. md.output.coordinates['positions'] to access final atom positions.
        """
        return None


def mdrun(input=None) -> MDRun:
    """MD simulation operation.

    Arguments:
        input : valid simulation input

    Returns:
        runnable operation to perform the specified simulation

    The returned object has a `run()` method to launch the simulation.
    Otherwise, this operation does not yet support the gmxapi data flow model.

    `input` may be a TPR file name.

    Note:
        New function names will be appearing to handle tasks that are separate

        "simulate" is plausibly a dispatcher or base class for various tasks
        dispatched by mdrun. Specific work factories are likely "minimize,"
        "test_particle_insertion," "legacy_simulation" (do_md), or "simulation"
        composition (which may be leap-frog, vv, and other algorithms)
    """
    # If input looks like classic filename input, use the gmxapi 0.0.7
    # implementation.
    from gmxapi.simulation.context import get_context
    try:
        work = workflow.from_tpr(input)
        assert work is not None
        work_context = get_context(work)
        assert work_context is not None
    except exceptions.UsageError:
        # Input was not gmxapi 0.0.7 input.
        work = None
        work_context = None

    # TODO: Inspect input to determine context.
    kwargs = {}
    handle = mdrun_dispatcher(context=work_context, input=work, label=None, **kwargs)

    return handle


class ReadTpr(OperationDetailsBase):
    """Produce simulation input data sources from a TPR file.

    Outputs:
        parameters : simulation parameters
        topology : molecular force field data
        coordinates : phase space coordinates
        simulation_state : integrator state and other checkpoint data

    The output is suitable as input for operations such as modify_input or mdrun.
    """
    _output = OutputCollectionDescription(parameters=dict)
    _input = InputCollectionDescription(
        [('filename', inspect.Parameter('filename',
                                        inspect.Parameter.POSITIONAL_OR_KEYWORD,
                                        annotation=str))])
    Resources = typing.Mapping

    def __call__(self, resources: Resources):
        # tpr_filename = resources.
        pass

    @classmethod
    def signature(cls) -> InputCollectionDescription:
        return cls._input

    @classmethod
    def resource_director(cls, *, input, output) -> Resources:
        pass

    def output_description(self) -> OutputCollectionDescription:
        return self._output

    def publishing_data_proxy(self, *, instance, client_id) -> DataProxyBase:
        pass

    def output_data_proxy(self, instance) -> DataProxyBase:
        pass

    def make_datastore(self, ensemble_width: int) -> typing.Mapping[str, OutputData]:
        description = ResultDescription(dtype=dict, width=ensemble_width)
        empty_output = OutputData(name='parameters', description=description)
        return {'parameters': empty_output}

    @classmethod
    def make_uid(cls, input) -> str:
        # TODO: input.fingerprint()
        # TODO: handle edge width
        filename = input.source_collection['filename']
        return 'read_tpr_' + filename


def read_tpr(filename=None):
    context = current_context()
    return ReadTpr.operation_director(filename=filename, context=context)

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
