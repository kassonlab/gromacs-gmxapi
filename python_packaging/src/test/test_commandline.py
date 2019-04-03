#!/usr/bin/env python
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

"""Unit tests for the gmxapi.commandline module."""

import unittest

from gmxapi import commandline
from gmxapi import util


class SimpleCliTestCase(unittest.TestCase):
    """Test creation and execution of the basic cli() command line wrapper."""

    def test_true(self):
        command = util.which('true')
        operation = commandline.cli(command=[command], shell=False)

        # Note: 'stdout' and 'stderr' not mapped.
        # Note: getitem not implemented.
        # assert not 'stdout' in operation.output
        # assert not 'stderr' in operation.output
        assert not hasattr(operation.output, 'stdout')
        assert not hasattr(operation.output, 'stderr')

        # Check for the attributes that we _do_ expect.
        assert hasattr(operation.output, 'erroroutput')
        assert hasattr(operation.output, 'returncode')

        operation.run()
        # assert operation.output.returncode.result() == 0
        assert operation.output.returncode == 0

    def test_false_explicit(self):
        command = util.which('false')
        operation = commandline.cli(command=[command], shell=False)
        # Explicitly run the operation.
        operation.run()
        assert operation.output.returncode == 1

    def test_false_implicit(self):
        command = util.which('false')
        operation = commandline.cli(command=[command], shell=False)
        operation.run()
        assert operation.output.returncode == 1

    def test_echo(self):
        # TODO: (FR5+) do we want to pipeline or checkpoint stdout somehow?
        operation = commandline.cli(command=[util.which('echo'), 'hi', 'there'])
        operation.run()
        assert operation.output.returncode == 0


if __name__ == '__main__':
    unittest.main()
