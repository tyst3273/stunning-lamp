"""
note: still need to add logic for it to check if the current phase space coordinate is the solution 
or not. actually, need to check if any of the moves for current coord are. if so, move there.
"""

import numpy as np

# --------------------------------------------------------------------------------------------------

class config:

    def __init__(self,coords,ensemble):

        """
        hold block coords for this particular phase space coordinate. find neibhoring 
        configurations and check if this config is connected to the solution
        """

        # this particular configs phase space coordinates
        self.coords = coords
        self.num_blocks = ensemble.num_blocks
        self.blocks = ensemble.blocks
        self.shapes = ensemble.shapes
        self.puzzle_shape = ensemble.puzzle_shape
        self.solution_block = ensemble.solution_block

        # the one to check against the solution
        self.masked = self.coords*(self.coords == self.solution_block).astype(int)- \
                            (self.coords != self.solution_block).astype(int)

        # check if we can move, check if is the solution 
        self._check_moves(ensemble.solution,ensemble.history)

    # ---------------------------------------------------------------------------------------------

    def _check_moves(self,solution,history):
        
        """
        find all neighboring configurations. not that since the coord we are moving from is 
        in the history, we wont go back to it
        """

        self.moves = []
        self.can_move = False
        self.solved = False
        self.solution = None

        # loop over the blocks and see if they can be moved
        for block in self.blocks:

            # spaces occupied by this block
            inds = np.argwhere(self.coords == block)
            rows = inds[:,0].reshape(self.shapes[block])
            cols = inds[:,1].reshape(self.shapes[block])

            # check moves for this block
            self._check_move('down',inds,rows,cols,block,history,solution)
            self._check_move('up',inds,rows,cols,block,history,solution)
            self._check_move('right',inds,rows,cols,block,history,solution)
            self._check_move('left',inds,rows,cols,block,history,solution)

    # ------------------------------------------------------------------------------------------

    def _check_move(self,direction,inds,rows,cols,block,history,solution):

        """
        check if the any moves are possible
        """

        if direction == 'down':
            shift = [1,0]
            bound = rows.max()
            shape = self.puzzle_shape[0]-1
        elif direction == 'up':
            shift = [-1,0]
            bound = rows.min() 
            shape = 0
        elif direction == 'right':
            shift = [0,1]
            bound = cols.max()
            shape = self.puzzle_shape[1]-1
        elif direction == 'left':
            shift = [0,-1]
            bound = cols.min()
            shape = 0
    
        if bound == shape: # check if at the top, dont move if so
            pass
        elif np.unique(self.coords[rows+shift[0],cols+shift[1]]).size != 1: # check if not empty
            pass
        elif np.unique(self.coords[rows+shift[0],cols+shift[1]])[0] != -1: # check if not empty  
            pass

        # otherwise do move
        else:
            move = np.copy(self.coords)
            move[rows,cols] = -1
            move[rows+shift[0],cols+shift[1]] = block

            # test if this configuration has been visited before
            flag = True
            for hist in history:
                test = hist-move
                if np.all(test): # if so, make note and quit checking
                    flag = False
                    break
            if flag: 
                self.moves.append(move)
                self.can_move = True

            # now test if this move is the solution
            test = solution-self.masked
            self.solved = np.all(test)
            if self.solved:
                self.solution = move

    # ----------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------

class ensemble:

    def __init__(self,init_coords=[[    0,  1,  1,  2],
                                   [    0,  1,  1,  2],
                                   [    3,  4,  4,  5],
                                   [    3,  6,  7,  5],
                                   [    8, -1, -1,  9]],
                         solution=[[   -1, -1, -1, -1],
                                   [   -1, -1, -1, -1],
                                   [   -1, -1, -1, -1],
                                   [   -1,  1,  1, -1],
                                   [   -1,  1,  1, -1]]):

        """
        ensemble of phase space coords visited by the solver. init_coords are the initial 
        coords of the blocks, solution is the configuration that solves it
        """
    
        # convert to arrays and error check shapes
        self.solution = np.array(solution)
        self.coords = np.array(init_coords)
        self.puzzle_shape = self.coords.shape
        if self.puzzle_shape != self.solution.shape:
            exit('fuck! shapes dont match')

        # print matrix
        print('\ninit:\n',self.coords,'\n\n')
        print('\nsolution:\n',self.solution,'\n\n')

        # record history of phase space coords
        self.history = []

        # get number and shapes of blocks
        self.num_blocks, self.blocks, self.shapes = self.get_blocks(self.coords)

        # check solution makes sense
        _n, self.solution_block, _s = self.get_blocks(self.solution)
        if _n != 1:
            exit('fuck! no solution given')
        self.solution_block = self.solution_block[0]
        _s = _s[0]
        ind = np.argwhere(self.blocks == self.solution_block).flatten()
        if ind.size == 0:
            exit('fuck! block doesnt exist')
        if _s != self.shapes[ind[0]]:
            exit('fuck! shape of solution is wrong')

    # ---------------------------------------------------------------------------------------------

    def solve(self):

        """
        use recursion and backtracking to solve the thing
        """

        # set initial state
        self.is_solved = False
        self.step = 0
        if not self.configs[0].can_move:
            exit('fuck! cant move')

        # this is recursive (i.e. solves itself)
        self.last = np.copy(self.coords)
        self.history.append(self.coords)
        self.config = config(self.coords,self)
        self.move()

        # check if it was solved
        if not self.is_solved:
            exit('fuck! it cant be solved')
        else:
            print(f'solved it willy nilly in {self.num_moves} moves')
            exit()

    # ---------------------------------------------------------------------------------------------

    def move(self):

        """
        use the current coords to find the next ones. DO NOT move DOWN into a coord we have already 
        visited. instead, if no moves are possible, move back up to the previous coord and visit its 
        next neighbor
        """

        print(f'\n\tnow on step {self.step+1}\n',self.config.coords)

        # check if this coord is connected to the solution
        if self.config.solved:
            print('\n\tsolved it!\n',self.config.solution)
            exit()

        # check if coords we just came UP into can move into new coord
        if len(self.configs[self.step].moves) == 0:
            return
        else:

            # the move to check
            move = self.configs[self.step].moves[0]

            # remove this 'move' from the possible ones since it is being added to the history
            self.configs[self.step].moves.pop(0)

            # test if this config has been visited before
            flag = True
            for hist in self.history:
                print(hist)
                test = hist-move
                if np.all(test):
                    flag = False
                    break

            # if not, create an new a phase space coord and move into it
            if flag:
                self.history.append(move)
                self.configs.append(config(move,self))
                self.step = self.step+1
                self.move()
        
    # ---------------------------------------------------------------------------------------------

    def get_blocks(self,coords):

        """
        get number of blocks, their shapes, etc.
        """

        # dont include < 0, it is empty space
        mask = (coords < 0).astype(int)
        coords = coords-mask*(coords+1)
        blocks = np.unique(coords)[1:]
        num_blocks = blocks.size

        shapes = []
        for block in blocks:
            inds = np.argwhere(coords == block)
            shapes.append([np.unique(inds[:,0]).size,np.unique(inds[:,1]).size])

        return num_blocks, blocks, shapes

    # -------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    ens = ensemble()
    ens.solve()
    



