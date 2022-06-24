#!/bin/bash

# '''
# Copyright 2011 Marco Mina. All rights reserved.
# 
# This file is part of fastSemSim
# 
# fastSemSim is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# fastSemSim is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with fastSemSim.  If not, see <http://www.gnu.org/licenses/>.
# '''

# set LOCAL to 0 to run the installed version of the GUI
# LOCAL=1

DIR_CODE=`dirname $0`
# echo $DIR_CODE

# Running the local version
python -c "from fastsemsim import fastsemsim_cmdline; fastsemsim_cmdline.start()" "$@"
