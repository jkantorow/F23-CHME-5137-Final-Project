"""
A module for the functions involved in building our project for simulating
phenolic resin curing in the presence of coal catalyzed by heat.

We are thinking of using a polymer ising model with kinetic Monte Carlo to
simulate the curing process. We will use a lattice to represent the polymer
components and the coal in an inital state.
"""

import numpy as np
import random as rand

def initial_random_matrix(n, ratio):

    """
    Generate an n x n matrix of randomly placed molecules
    with a user-defined distrubution ratio.
    We will define "molecules" as a phenol, a coal molecule, 
    and, although not a physical object, an empty space.
    All of these "molecules" will be of comparable size for
    simplicity.
    """

    state = np.zeros((n, n))

    # We can define:

    mol_dict = {
        0 : "void",
        1 : "phenol",
        2 : "coal"
    }
    
    # Ratio of empty space to phenol to coal in the form
    # (x:y:z) where x, y, and z are integers:

    void_ratio = ratio[0]
    phenol_ratio = ratio[1]
    coal_ratio = ratio[2]

    # Randomize each cell of the state matrix according to the 
    # user-defined ratio:

    for i in range(n):
        for j in range(n):
            r = rand.randint(1, void_ratio + phenol_ratio + coal_ratio)
            if r <= void_ratio:
                state[i][j] = 0
            elif r <= void_ratio + phenol_ratio:
                state[i][j] = 1
            else:
                state[i][j] = 2
    
    return state


# Define and get the rates of each event occuring the system based on the
# input parameters:

def get_rates(grid, mol_type, pos):

    """
    We don't necessarily know what the rates should be for each of these
    events; for now we use a place holder to indicate that we need to 
    update this.

    We need to determine what the input parameters should be
    to describe each molecule of the system (we can change these as needed):

    grid     : numpy array - the current state of the system.

    mol_type : float - the type of molecule in the system (1, or 2);
             : all information about the molecule can be stored in this
             : function and used to determine the rates of each event.

    pos      : tuple - the position of the molecule in the matrix; this
             : is used to find the identification of the peripheral molecules
             : in the immediate vicinity of the current molecule.

             molecule i,j and its immediate periphery molecules can be
             thought of as follows where i=rows and j=columns:

             (i-1,j-1) (i-1,j  ) (i-1,j+1)

             (i  ,j-1) (i  ,j  ) (i  ,j+1)

             (i+1,j-1) (i+1,j  ) (i+1,j+1)

    """

    # Get identities and positions of each peripheral molecule:
    # by using a 3x3 mask:

    x, y = pos
    periphery = []

    for i in range(x-1, x+2):
        for j in range(y-1, y+2):
            periphery.append(grid[i][j])

    periphery = np.array(periphery).reshape(3, 3)

    # Calculate the rate of each event occuring based
    # on a huge if else tree:

    # If mol_type is a phenol:

    if mol_type == 1:

        for i in range(3):
            for j in range(3):

                # if (i, j) is the central molecule, skip it:

                if i == x and j == y:
                    continue

                elif periphery[i, j] == 1:
                    PLACEHOLDER # NEED EDIT

    PLACEHOLDER = None

    # Reaction events:

    # KMC rate of phenol reacting with another phenol:
    # I think we can determine this rate mostly using the
    # structure of the molecules interacting and the
    # reaction rate constant (theoretically).

    rate_res_cures_alone = PLACEHOLDER

    # KMC rate of phenol reacting with coal:

    rate_res_cures_with_coal = PLACEHOLDER

    # KMC rate of no reaction occuring:

    rate_no_rxn = PLACEHOLDER

    # Translocation events (could nix this if
    # the system gets too complicated):

    # Getting rid of the swapping problem could
    # make this much simpler, also Adam suggests
    # compiling all of these rates into one "
    # movement" rate.

    # To do this let's define a singular rate
    # of a movement event occuring:

    rate_mol_moves = PLACEHOLDER

    rate_mol_moves_NW = PLACEHOLDER
    rate_mol_moves_N = PLACEHOLDER
    rate_mol_moves_NE = PLACEHOLDER
    rate_mol_moves_E = PLACEHOLDER
    rate_mol_moves_SE = PLACEHOLDER
    rate_mol_moves_S = PLACEHOLDER
    rate_mol_moves_SW = PLACEHOLDER
    rate_mol_moves_W = PLACEHOLDER


    # We need the output of the crosslinking degree
    # depending on the molecular weight and the
    # involved functional group of the molecule.
    # We can keep track of every time a crosslinking
    # event occurs and then represent this value as a ratio
    # of the total number of molecules in the system.

    return (rate_res_cures_alone,
            rate_res_cures_with_coal,
            rate_no_rxn,
            rate_mol_moves_NW,
            rate_mol_moves_N,
            rate_mol_moves_NE,
            rate_mol_moves_E,
            rate_mol_moves_SE,
            rate_mol_moves_S,
            rate_mol_moves_SW,
            rate_mol_moves_W)

