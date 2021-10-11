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
import numpy as np

# custom modules
import modpack.printing as prt
import modpack.error as err
import modpack.warning as warn

# --------------------------------------------------------------------------------------------------

class params:

    """
    this class is designed to get and hold all input variables from the user. it should really
    only be used to define input parameters, i.e. everything that would go in an 'input file' 
    for a DFT code. should provide as much error checking as possible in here and return, as much
    as possible, the final formatted version of the variable expected downstream in the code. 
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self):

        """
        define 'status' of variable and all defaults here
        """

        # variables required by the user
        self.required =  ['lat_vecs',
                          'lat_scale',
                          'atom_pos']

        # optional variables
        self.optional =  ['tasks',
                          'pw_cutoff',
                          'lo_cutoff',
                          'verbosity',
                          'xc_type',
                          'kpoints',
                          'wf_in',
                          'wf_out',
                          'den_in',
                          'den_out']

        # variables not ready for production use
        self.developer = ['xc_type',
                          'kpoints',
                          'wf_in',
                          'wf_out',
                          'den_in',
                          'den_out']

        # a set of all allowed variables
        self.allowed = self.optional+self.required

        # default iterable variables
        self.tasks =['scf']
        self.pw_cutoff = [50]
        self.lo_cutoff = [100]
        self.kpoints = ['from_file']
        self.wf_in = ['atomic']
        self.wf_out = ['wf_1']
        self.den_in = ['from_wf']
        self.den_out = ['den_1']

        # default global variables
        self.verbosity = 1
        self.xc_type = 'simple'

    # ----------------------------------------------------------------------------------------------

    def parse_params(self,input_file,parse_method='text'):

        """
        method to parse variables from a text file. right now, only text files are supported,
        but i plan to make hdf5 and maybe yaml or something supported
        """

        # see if input file exists
        if not os.path.exists(input_file):
            err.file_error(input_file)

        # see if input file works
        self.input_file = input_file
        try:
            with open(self.input_file,'r') as f_in:
                self.input_lines = f_in.readlines()
                self.input_text = f_in.read()
        except:
            err.file_error(self.input_file)

        # get text from input files and process it
        self._preprocess_text()
        self._check_vars()

        # get the variables from the variable dictionary 
        self._get_all_vars()

        # now (possibly) print the variables to the screen
        if self.verbosity > 0:
            self._echo_variables()

    # ----------------------------------------------------------------------------------------------

    def _preprocess_text(self):

        """
        get text from file, strip comment and blank lines, create dictionary of variables 
        called 'var_dict'. this is the essential result of this method and contains all keywords
        that are allowed in the input with the arguments given to them as key-value pairs.
        """

        keyword_inds = [] # position of keywords
        keywords = [] # the keywords found in the file 
        ind = 0

        trimmed_lines = []
        for line in self.input_lines:   

            # skip blank and comment lines
            if line.strip().startswith('#') or not len(line.strip().split()):
                continue

            else:
                line = line.split('#')[0].strip() # strip comments
                trimmed_lines.append(line)

                # check if its a keyword
                if line.endswith(':'):
                    keyword_inds.append(ind)
                    keywords.append(line.strip(':'))

                ind = ind+1

        self.input_lines = trimmed_lines
        num_kw = len(keywords)

        # now put it all into the variable dictionary
        self.var_dict = {}
        for ii in range(num_kw-1):
            self.var_dict[keywords[ii]] = self.input_lines[keyword_inds[ii]+1:keyword_inds[ii+1]]
        self.var_dict[keywords[-1]] = self.input_lines[keyword_inds[-1]+1:]

        # now get the 'unique' list of kw from the input file
        self.keywords = self.var_dict.keys()

        # get full list of keywords, wrong and everything
        self.raw_keywords = keywords 

    # ----------------------------------------------------------------------------------------------

    def _check_vars(self):

        """
        if not all required variables or any unkown variables are given, or if any variables
        are duplicated crash. print a warned about all developer variables
        """

        unknown = []
        checked = []
        duplicates = []
        dev = []
        empty = []

        # check the the variables
        for kw in self.raw_keywords:    
            if kw not in self.allowed:
                unknown.append(kw)
            elif kw in self.required:
                ind = self.required.index(kw)
                self.required.pop(ind)
            if kw in self.developer:
                dev.append(kw)
            if len(self.var_dict[kw]) == 0:
                empty.append(kw)
            if kw in checked:
                duplicates.append(kw)
            checked.append(kw)

        # crash if needed
        if len(self.required) != 0:
            err.missing_keyword_error(self.required)
        if len(unknown) != 0:
            err.unknown_keyword_error(unknown)
        if len(empty) != 0:
            err.empty_keyword_error(empty)
        if len(duplicates) != 0:
            err.duplicated_keyword_error(duplicates)

        # check if developer variables were used
        if len(dev) != 0:
            warn.developer_keyword_warning(dev)

    # ----------------------------------------------------------------------------------------------

    def _echo_variables(self):

        """
        echo the variables to the screen
        """

        # print info
        msg_type = 'INPUT VARIABLES'
        prt.print_stdout('echoing input variables to the screen',msg_type)

        # num_tasks
        msg = 'num_tasks:\n'
        msg = msg+f'  {self.num_tasks:<}'
        prt.print_stdout(msg,msg_type=None)

        # tasks
        msg = 'tasks:\n'
        for ii in range(self.num_tasks-1):
            msg = msg+f'  [{ii}] : {self.tasks[ii]:<}\n'
        msg = msg+f'  [{self.num_tasks-1}] : {self.tasks[-1]:<}'
        prt.print_stdout(msg,msg_type=None)
        
        # verbosity
        msg = 'verbosity:\n'
        msg = msg+f'  {self.verbosity:<}'
        prt.print_stdout(msg,msg_type=None)

        # num_atoms
        msg = 'num_atoms:\n'
        msg = msg+f'  {self.num_atoms:<}'
        prt.print_stdout(msg,msg_type=None)

        # atom_pos
        msg = 'atom_pos (crystal coords):\n'
        for ii in range(self.num_atoms):
            msg = msg+f'  ({ii}) : '
            for jj in range(3):
                msg = msg+f'{self.atom_pos[ii,jj]: <10.6f}'
            if ii == self.num_atoms-1:
                continue
            else:
                msg = msg+'\n'
        prt.print_stdout(msg,msg_type=None)

        # lat_vecs
        msg = 'lat_vecs (Angstrom):\n'
        for ii in range(3):
            msg = msg+'  '
            for jj in range(3):
                msg = msg+f'{self.lat_vecs[ii,jj]: <10.6f}'
            if ii == 2:
                continue
            else:
                msg = msg+'\n'
        prt.print_stdout(msg,msg_type=None)

        # lat_scale
        msg = 'lat_scale:\n'
        for ii in range(self.num_tasks-1):
            msg = msg+f'  [{ii}] : {self.lat_scale[ii]:<}\n'
        msg = msg+f'  [{self.num_tasks-1}] : {self.lat_scale[-1]:<}'
        prt.print_stdout(msg,msg_type=None)

        # pw_cutoff
        msg = 'pw_cutoff (eV):\n'
        for ii in range(self.num_tasks-1):
            msg = msg+f'  [{ii}] : {self.pw_cutoff[ii]:<}\n'
        msg = msg+f'  [{self.num_tasks-1}] : {self.pw_cutoff[-1]:<}'
        prt.print_stdout(msg,msg_type=None)

        # lo_cutoff
        msg = 'lo_cutoff (eV):\n'
        for ii in range(self.num_tasks-1):
            msg = msg+f'  [{ii}] : {self.lo_cutoff[ii]:<}\n'
        msg = msg+f'  [{self.num_tasks-1}] : {self.lo_cutoff[-1]:<}'
        prt.print_stdout(msg,msg_type=None)

    # ----------------------------------------------------------------------------------------------

    def _get_all_vars(self):

        """
        get the variables from the input file. should provide as much error checking as possible 
        in here and return, as much as possible, the final formatted version of the variable 
        expected downstream in the code. 

        the variable num_tasks is created based on the number of tasks to be done.
        variables that are 'iterable' must have either 1 argument or num_tasks arguments. if only
        1, it is copied num_tasks times. 'global' variables are not copied and the same values 
        are used for all tasks (unless the code modifies them!)

        if new variables added that depend on others, be careful to add them in the logically 
        consistent place
        """

        self._get_tasks()
        self._get_verbosity()
        self._get_lat_scale()
        self._get_lat_vecs()
        self._get_pw_cutoff()
        self._get_lo_cutoff()
        self._get_atom_pos()

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # get the variables and error check them

    def _get_tasks(self):

        """
        tasks variable. see documentation in 'meta' (dictionary of metadata) below
        """

        meta = {'name':
                    'tasks',
                'intent':
                    'iterable',
                'input format':
                    'list of strings',
                'description':
                    'list of tasks to be done by the code. the code loops over them, \n '
                    'doing each in turn.',
                'allowed args':
                    '- scf \n '
                    '  do scf procedure to convege the density to user specified value \n '
                    '- nscf \n '
                    '  do non-scf procedure to get wavefunctions. if not using exact \n '
                    '  diagonalization, the convergence criterion should be the wave \n '
                    '  function residuals.'}

        # get from var_dict
        self.tasks = self._parse_string_list('tasks',self.tasks,meta=meta)

        # check that they are allowed:
        for task in self.tasks:
            if task not in ['scf','nscf']:
                err.args_wrong_error('tasks',self.tasks,meta)

        # the number of tasks to be done
        self.num_tasks = len(self.tasks)

    # ----------------------------------------------------------------------------------------------

    def _get_verbosity(self):

        """
        verbosity variable. see documentation in 'meta' (dictionary of metadata) below
        """

        meta = {'name':
                    'verbosity',
                'intent':
                    'global',
                'input format':
                    'single integer',
                'description':
                    'controls how much info is printed to the screen',
                'allowed args':
                    '- 0: print the minimum amount \n'
                    ' - 1: the usual amount of printing \n'
                    ' - 2: print everything \n'}

        # get from var_dict
        self.verbosity = self._parse_int_list('verbosity',self.verbosity, \
                    meta,num_rows=[1],num_col=[1])[0][0] # single int
    
        # check that they are allowed:
        if self.verbosity not in [0,1,2]:
            err.args_wrong_error('verbosity',[self.verbosity],meta)

    # ----------------------------------------------------------------------------------------------

    def _get_lat_scale(self):

        """
        lat_scale. see documentation in 'meta' (dictionary of metadata) below
        """

        meta = {'name':
                    'lat_scale',
                'intent':
                    'iterable',
                'input format':
                    'list of floats',
                'description':
                    'scaling factor for lattice vectors',
                'allowed args':
                    'a float that is > 0'}

        # get from var_dict
        self.lat_scale = self._parse_float_list('lat_scale',default=None, \
                meta=meta,num_rows=[1,self.num_tasks],num_col=[1])
        self.lat_scale = np.array(self.lat_scale).flatten()

        # copy to match number of tasks
        if self.lat_scale.size == 1:
            self.lat_scale = np.tile(self.lat_scale,reps=[self.num_tasks])

        # check that they are allowed:
        for arg in self.lat_scale:
            if not arg > 0:
                err.args_wrong_error('lat_scale',self.lat_scale,meta)

    # ----------------------------------------------------------------------------------------------

    def _get_lat_vecs(self):

        """
        lat_vecs. see documentation in 'meta' (dictionary of metadata) below
        """

        meta = {'name':
                    'lat_vecs',
                'intent':
                    'global',
                'input format':
                    '3x3 list of floats',
                'description':
                    'lattice vectors in angstroms (row vectors)'}

        # get from var_dict
        self.lat_vecs = self._parse_float_list('lat_vecs',default=None, \
                meta=meta,num_rows=[3],num_col=[3])
        self.lat_vecs = np.array(self.lat_vecs)

    # ----------------------------------------------------------------------------------------------

    def _get_pw_cutoff(self):

        """
        pw_cutoff. see documentation in 'meta' (dictionary of metadata) below
        """

        meta = {'name':
                    'pw_cutoff',
                'intent':
                    'iterable',
                'input format':
                    'list of floats',
                'description':
                    'energy cutoff for the expansion of the pw part of the basis set (in eV)',
                'allowed args':
                    'a list of floats, all > 0'}

        # get from var_dict
        self.pw_cutoff = self._parse_float_list('pw_cutoff',default=self.pw_cutoff, \
                meta=meta,num_rows=[1,self.num_tasks],num_col=[1])
        self.pw_cutoff = np.array(self.pw_cutoff).flatten()

        # copy to match number of tasks
        if self.pw_cutoff.size == 1:
            self.pw_cutoff = np.tile(self.pw_cutoff,reps=[self.num_tasks])

        # check that they are allowed:
        for arg in self.pw_cutoff:
            if not arg > 0:
                err.args_wrong_error('pw_cutoff',self.pw_cutoff,meta)

    # ----------------------------------------------------------------------------------------------

    def _get_lo_cutoff(self):

        """
        lo_cutoff. see documentation in 'meta' (dictionary of metadata) below
        """

        meta = {'name':
                    'lo_cutoff',
                'intent':
                    'iterable',
                'input format':
                    'list of floats',
                'description':
                    'energy cutoff for the expansion of the lo part of the basis set (in eV) \n',
                    ' please note that the lo part is expanded in pw\'s too.'
                'allowed args':
                    'a list of floats, all > 0'}

        # get from var_dict
        self.lo_cutoff = self._parse_float_list('lo_cutoff',default=self.lo_cutoff, \
                meta=meta,num_rows=[1,self.num_tasks],num_col=[1])
        self.lo_cutoff = np.array(self.lo_cutoff).flatten()

        # copy to match number of tasks
        if self.lo_cutoff.size == 1:
            self.lo_cutoff = np.tile(self.lo_cutoff,reps=[self.num_tasks])

        # check that they are allowed:
        for arg in self.lo_cutoff:
            if not arg > 0:
                err.args_wrong_error('lo_cutoff',self.lo_cutoff,meta)

    # ---------------------------------------------------------------------------------------------

    def _get_atom_pos(self):

        """
        atom_pos. see documentation in 'meta' (dictionary of metadata) below
        """

        meta = {'name':
                    'atom_pos',
                'intent':
                    'global',
                'input format':
                    'natom*3 list of floats',
                'description':
                    'atom positions in crystal coordinates'}

        # get from var_dict
        self.atom_pos = self._parse_float_list('atom_pos',default=None, \
                meta=meta,num_rows=None,num_col=[3])
        self.atom_pos = np.array(self.atom_pos)

        # get the number of atoms
        self.num_atoms = self.atom_pos.shape[0]

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # methods to get the data from the data structure

    def _parse_string_list(self,keyword,default=None,meta=None,num_rows=None):

        """
        return a list of strings. only 1 column is expected
        """

        # see if keyword is present
        if keyword not in self.keywords:
            return default

        # get it from the dictionary
        args = self.var_dict[keyword]
        num_args = len(args)

        # strip any quotation marks
        args = [arg.strip().strip('\'"') for arg in args]

        # check that the number of args is expected
        if num_rows != None:
            if num_args not in num_rows:
                err.args_wrong_error(keyword,args,meta)

        # only one column per row
        for arg in args:
            if len(arg.split()) != 1:
                err.args_wrong_error(keyword,args,meta)
        
        return args

    # ----------------------------------------------------------------------------------------------

    def _parse_int_list(self,keyword,default=None,meta=None,num_rows=None,num_col=None):

        """
        return a list of ints
        """

        # see if keyword is present
        if keyword not in self.keywords:
            return default

        # get it from the dictionary
        args = self.var_dict[keyword]
        num_args = len(args)

        # convert to ints
        try:
            for ii in range(num_args):
                args[ii] = [int(arg.strip()) for arg in args[ii].split()]
        except:
            err.args_wrong_error(keyword,args,meta)
            
        # check that the number of args is expected
        if num_rows != None:
            if num_args not in num_rows:
                err.args_wrong_error(keyword,args,meta)

        # check that number of columns is as expected
        if num_col != None:
            for arg in args:
                if len(arg) not in num_col:
                    err.args_wrong_error(keyword,args,meta)

        return args

    # ----------------------------------------------------------------------------------------------

    def _parse_float_list(self,keyword,default=None,meta=None,num_rows=None,num_col=None):

        """
        return a list of floats
        """

        # see if keyword is present
        if keyword not in self.keywords:
            return default

        # get it from the dictionary
        args = self.var_dict[keyword]
        num_args = len(args)

        # convert to floats
        try:
            for ii in range(num_args):
                args[ii] = [float(arg) for arg in args[ii].split()]
        except:
            err.args_wrong_error(keyword,args,meta)

        # check that the number of args is expected
        if num_rows != None:
            if num_args not in num_rows:
                err.args_wrong_error(keyword,args,meta)

        # check that number of columns is as expected
        if num_col != None:
            for arg in args:
                if len(arg) not in num_col:
                    err.args_wrong_error(keyword,args,meta)

        return args

# --------------------------------------------------------------------------------------------------










