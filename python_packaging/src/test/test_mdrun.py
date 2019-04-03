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

"""Test gromacs.mdrun operation.

Factory produces deferred execution operation.

Factory expects input for conformation, topology, simulation parameters, simulation state.

TODO: Factory accepts additional keyword input to indicate binding
 to the "potential" interface.
"""

import json
import logging
import os
import pytest
import shutil
import tempfile
import unittest
import warnings

# Configure the `logging` module before and non-built-in packages start to use it.
logging.getLogger().setLevel(logging.DEBUG)
# create console handler
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
# create formatter and add it to the handler
formatter = logging.Formatter('%(asctime)s:%(name)s:%(levelname)s: %(message)s')
ch.setFormatter(formatter)
# add the handlers to the logger
logging.getLogger().addHandler(ch)

import gmxapi as gmx

from pytesthelpers import withmpi_only

@pytest.mark.usefixtures('cleandir')
def test_run_from_tpr(spc216):
    assert os.path.exists(spc216)

    # sim_input = gmx.read_tpr(self.tprfilename)
    # md = gmx.mdrun(sim_input)
    # md.run()

    md = gmx.mdrun(spc216)
    md.run()


# @pytest.mark.usefixtures("cleandir")
# class BindingsTestCase(unittest.TestCase):
#     def setUp(self):
#         self.tpr_filename
#     def test_APIObjectsFromTpr(self):
#         apisystem = core.from_tpr(tpr_filename)
#         assert isinstance(apisystem, core.MDSystem)
#         context = core.Context()
#         mdargs = core.MDArgs()
#         mdargs.set({'threads_per_rank': 1})
#         context.setMDArgs(mdargs)
#         assert hasattr(apisystem, 'launch')
#         session = apisystem.launch(context)
#         assert hasattr(session, 'run')
#         session.run()
#
# @pytest.mark.usefixtures("cleandir")
# @pytest.mark.usefixtures("caplog")
# def test_simpleSimulation(caplog):
#     """Load a work specification with a single TPR file and run."""
#     # use case 1: simple high-level
#     md = gmx.workflow.from_tpr(tpr_filename, threads_per_rank=1)
#     gmx.run(md)
#
# @pytest.mark.usefixtures("cleandir")
# @pytest.mark.filterwarnings("ignore:Using or importing the ABCs from 'collections'")
# @pytest.mark.usefixtures("caplog")
# def test_idempotence1(caplog):
#     """Confirm that a work graph can be run repeatedly, even after completed.
#
#     Use gmx.run and avoid extra references held by user code.
#     """
#     md = gmx.workflow.from_tpr(tpr_filename, threads_per_rank=1)
#     gmx.run(md)
#     gmx.run(md)
#     gmx.run(md)
#     md = gmx.workflow.from_tpr(tpr_filename, threads_per_rank=1)
#     gmx.run(md)
#     md = gmx.workflow.from_tpr(tpr_filename, threads_per_rank=1)
#     gmx.run(md)
#
# @pytest.mark.usefixtures("cleandir")
# @pytest.mark.filterwarnings("ignore:Using or importing the ABCs from 'collections'")
# @pytest.mark.usefixtures("caplog")
# def test_idempotence2(caplog):
#     """Confirm that a work graph can be run repeatedly, even after completed.
#
#     Interact with Context more directly.
#     Check that more unpredictable references held by user are still safe.
#     """
#     md = gmx.workflow.from_tpr(tpr_filename, threads_per_rank=1)
#     with gmx.get_context(md) as session:
#         session.run()
#
#     context = gmx.get_context(md)
#     with context as session:
#         session.run()
#
#     context = gmx.context.Context()
#     context.work = md
#     with context as session:
#         session.run()
#
# @pytest.mark.usefixtures("cleandir")
# def test_modifiedInput(caplog):
#     """Load a work specification with a single TPR file and updated params."""
#     md = gmx.workflow.from_tpr(tpr_filename, threads_per_rank=1, end_time='0.02')
#     context = gmx.get_context(md)
#     with context as session:
#         session.run()
#     md = gmx.workflow.from_tpr(tpr_filename, threads_per_rank=1, end_time='0.03')
#     context.work = md
#     with context as session:
#         session.run()
#     md = gmx.workflow.from_tpr(tpr_filename, threads_per_rank=1, end_time='0.04')
#     gmx.run(md)
#
# @pytest.mark.usefixtures("cleandir")
# @pytest.mark.usefixtures("caplog")
# @withmpi_only
# def test_plugin_no_ensemble(caplog):
#     # Test attachment of external code
#     md = gmx.workflow.from_tpr(tpr_filename, threads_per_rank=1)
#
#     # Create a WorkElement for the potential
#     #potential = gmx.core.TestModule()
#     potential_element = gmx.workflow.WorkElement(namespace="testing", operation="create_test")
#     potential_element.name = "test_module"
#     before = md.workspec.elements[md.name]
#     md.add_dependency(potential_element)
#     assert potential_element.name in md.workspec.elements
#     assert potential_element.workspec is md.workspec
#     after = md.workspec.elements[md.name]
#     assert not before is after
#
#     # Workaround for https://github.com/kassonlab/gmxapi/issues/42
#     # We can't add an operation to a context that doesn't exist yet, but we can't
#     # add a work graph with an operation that is not defined in a context.
#     context = gmx.get_context()
#     context.add_operation(potential_element.namespace, potential_element.operation, my_plugin)
#     context.work = md
#
#     with warnings.catch_warnings():
#         # Swallow warning about wide MPI context
#         warnings.simplefilter("ignore")
#         with context as session:
#             if context.rank == 0:
#                 print(context.work)
#             session.run()
#
#
# @pytest.mark.usefixtures("cleandir")
# @pytest.mark.usefixtures("caplog")
# @withmpi_only
# def test_plugin_with_ensemble(caplog):
#     # Test in ensemble.
#     md = gmx.workflow.from_tpr([tpr_filename, tpr_filename], threads_per_rank=1)
#
#     # Create a WorkElement for the potential
#     #potential = gmx.core.TestModule()
#     potential_element = gmx.workflow.WorkElement(namespace="testing", operation="create_test")
#     potential_element.name = "test_module"
#     before = md.workspec.elements[md.name]
#     md.add_dependency(potential_element)
#     assert potential_element.name in md.workspec.elements
#     assert potential_element.workspec is md.workspec
#     after = md.workspec.elements[md.name]
#     assert not before is after
#
#     # Workaround for https://github.com/kassonlab/gmxapi/issues/42
#     # We can't add an operation to a context that doesn't exist yet, but we can't
#     # add a work graph with an operation that is not defined in a context.
#     context = gmx.get_context()
#     context.add_operation(potential_element.namespace, potential_element.operation, my_plugin)
#     context.work = md
#
#     with warnings.catch_warnings():
#         # Swallow warning about wide MPI context
#         warnings.simplefilter("ignore")
#         with context as session:
#             if context.rank == 0:
#                 print(context.work)
#             session.run()


if __name__ == '__main__':
    unittest.main()
