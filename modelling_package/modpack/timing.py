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
from timeit import default_timer

# custom modules
import modpack.error as err
import modpack.printing as prt
import modpack.warning as wrn

# --------------------------------------------------------------------------------------------------

class _timer:

    """
    object for each 'timer'
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self,label,units,verbosity_vs):

        """
        holds the data for each timer
        """
        
        if units == 's':
            scale = 1
        elif units == 'm':
            scale = 1/60
        elif units == 'ms':
            scale = 1e3
        else:
            units = 's'
            scale = 1
            
        self.label = label
        self.start_time = 0
        self.end_time = 0
        self.lap_time = 0
        self.tot_time = 0
        self.units = units
        self.scale = scale
        self.verbosity_vs = verbosity_vs
        self.running = False
        self.num_calls = 0
    
    # ----------------------------------------------------------------------------------------------

    def _start(self):

        """
        start the timer and set its info
        """
        
        self.start_time = default_timer()
        self.num_calls = self.num_calls+1
        self.running = True

    # ----------------------------------------------------------------------------------------------

    def _stop(self):

        """
        stop the timer and update its info
        """

        self.end_time = default_timer()
        self.lap_time = (self.end_time-self.start_time)*self.scale
        self.tot_time = self.lap_time+self.tot_time
        self.running = False


# --------------------------------------------------------------------------------------------------

class timers:
    
    """
    holds 'timer' objects that are used to time portions of the code and print to the screen
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self,verbosity):

        """
        the dictionary maps timer objects to thier labels
        """

        self.timer_dict = {}
        self.verbosity = verbosity

    # ----------------------------------------------------------------------------------------------

    def start_timer(self,label,units='s',verbosity_vs=1):
        
        """
        if the timer doesnt exist, create and start it. otherwise, update the tot and lap time, 
        etc
        """

        # see if it exits
        if label not in self.timer_dict.keys():
            self.timer_dict[label] = _timer(label,units,verbosity_vs)

        # or if its running
        elif self.timer_dict[label].running:
            self.timer_dict[label]._stop()

        # now start it
        self.timer_dict[label]._start()


    # ----------------------------------------------------------------------------------------------

    def stop_timer(self,label):

        """
        stop the timer if it is running. if its not running, print a warning that it should be
        """

        # see if it exits
        if label not in self.timer_dict.keys():
            wrn.generic_warning(msg=f'the timer \'{label}\' doesnt exist!')
        else:
            
            # check if its running
            if self.timer_dict[label].running:
                self.timer_dict[label]._stop()
            else:
                wrn.generic_warning(msg=f'the timer \'{label}\' wasnt running!')

    # ----------------------------------------------------------------------------------------------

    def stop_all_timers(self):

        """
        stop all running timers
        """

        for label in self.timer_dict.keys():
            if self.timer_dict[label].running:
                self.timer_dict[label]._stop()

    # ----------------------------------------------------------------------------------------------

    def print_timing(self):

        """
        print all timing infor
        """

        prt.print_bar('%')
        msg = 'print timing information\n'
        prt.print_stdout(msg,msg_type='TIMING')
        msg = ' label               running  num_calls  lap_time  tot_time  units \n'
        msg = msg+'  -----               -------  ---------  --------  --------  -----\n' 
        for label in self.timer_dict.keys():

            # only print the timer if verbosity is high enough
            if self.verbosity < self.timer_dict[label].verbosity_vs:
                continue

            running = int(self.timer_dict[label].running)
            lap_time = self.timer_dict[label].lap_time
            tot_time = self.timer_dict[label].tot_time
            num_calls = self.timer_dict[label].num_calls
            units = self.timer_dict[label].units
            msg = msg+f'  {label:<20}     {running}    {num_calls:5g}    {lap_time:8.4f} ' \
                 f'  {tot_time:8.4f}   [{units}]\n'

        prt.print_stdout(msg)

# --------------------------------------------------------------------------------------------------












