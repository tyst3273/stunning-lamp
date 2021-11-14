# -*- coding: utf-8 -*-

#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   !                                                                                    !
#   ! Copyright (C) 2021 Tyler Sterling                                                  !
#   !  All rights reserved.                                                              !
#   !                                                                                    !
#   ! This file is part of modpack.                                                      !
#   !                                                                                    !
#   ! Redistribution and use in source and binary forms, with or without modification,   !
#   ! are permitted provided that the following conditions are met:                      !
#   !                                                                                    !
#   ! * Redistributions of source code must retain the above copyright notice, this list !
#   !   of conditions and the following disclaimer.                                      !
#   !                                                                                    !
#   ! * Redistributions in binary form must reproduce the above copyright notice, this   !
#   !   list of conditions and the following disclaimer in the documentation and/or      !
#   !   other materials provided with the distribution.                                  !
#   !                                                                                    !
#   ! * Neither the name of the modpack project nor the names of its contributors may    !
#   !   be used to endorse or promote products derived from this software without        !
#   !   specific prior written permission.                                               !
#   !                                                                                    !
#   ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    !
#   ! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      !
#   ! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. !
#   ! IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,   !
#   ! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT !
#   ! NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR !
#   ! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,  !
#   ! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) !
#   ! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE         !
#   ! POSSIBILITY OF SUCH DAMAGE.                                                        !
#   !                                                                                    !
#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# system modules
import os
import sys

# custom modules
import modpack.input_params as input_params
import modpack.timing as timing
import modpack.error as err
import modpack.printing as prt
import modpack.models.kohn_sham as ks

# --------------------------------------------------------------------------------------------------

# get input file (and command line args) from the user
if len(sys.argv) != 1:
    input_file = sys.argv[1]

# or use default
else:
    input_file = 'params.in'

# --------------------------------------------------------------------------------------------------

# get the input parameters from the file
params = input_params.params()
params.parse_params(input_file)

# start up the timers
timers = timing.timers(params.verbosity)

# --------------------------------------------------------------------------------------------------

# print the list of tasks to be done
msg = f'there are {params.num_tasks} tasks to do:\n'
for ii in range(params.num_tasks-1):
    msg = msg+f'  [{ii}] : \'{params.tasks[ii]}\'\n'
msg = msg+f'  [{params.num_tasks-1}] : \'{params.tasks[-1]} '
prt.print_stdout(msg)

# start outer loop timer
timers.start_timer('outer_loop')

# loop over all the tasks
for ii in range(params.num_tasks):

    task = params.tasks[ii]

    # time this task
    timers.start_timer(f'task_{ii}',units='ms')
    
    # print info
    prt.print_bar('-')
    msg = f'now doing task[{ii}]: {task}'
    prt.print_stdout(msg)

    # set up the kohn-sham model hamiltonian
    if task in ['scf','nscf']:
        ks_model = ks.ks_model(params,timers)

    # stop this timer
    timers.stop_timer(f'task_{ii}')

# --------------------------------------------------------------------------------------------------

# stop all timers and print final timing
timers.stop_all_timers()
timers.print_timing()

# --------------------------------------------------------------------------------------------------














