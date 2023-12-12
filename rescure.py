"""
A module for the functions involved in building our project for simulating
phenolic resin curing in the presence of coal catalyzed by heat.

We are thinking of using a polymer ising model with kinetic Monte Carlo to
simulate the curing process. We will use a lattice to represent the polymer
components and the coal in an inital state. This lattice will represent a 
cross-section of the system to capture the process of cross-linking.

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

    It's important to recognize that the 2D plane of this
    matrix represents the perpendicular cross section of
    multiple parallel resin polymers in close (angstroms)
    proximity and therefore the crosslinking of phenol /
    coal residues in the same perpendicular plane yet
    on different polymer chains.

    """

    state = np.zeros((n, n))

    # We can define:

    mol_dict = {
        0 : "void",
        1 : "phenol",
        2 : "coal"
    }

    # It is important to note that more numbers can be used
    # to represent crosslinked phenol molecules once the
    # simulation is run, although we do not start with any
    # crosslinked molecules in the initial state.
    
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

def get_rates_prototype(state, T, mol_type, pos):

    """
    We don't necessarily know what the rates should be for each of these
    events; for now we use a place holder to indicate that we need to 
    update this.

    We need to determine what the input parameters should be
    to describe each molecule of the system (we can change these as needed):

    state    : numpy array - the current state of the system in matrix
             : format.

    T        : float - the temperature of the molecule (same as the temperature
             : of the whole system in our case).

    mol_type : float - the type of molecule in the system (1, or 2);
             : all information about the molecule can be stored in this
             : function and used to determine the rates of each event.

    pos      : tuple - the position of the molecule in the matrix; this
             : is used to find the identification of the peripheral molecules
             : in the immediate vicinity of the current molecule.

             molecule i,j and its immediate peripheral molecules can be
             thought of as follows where i=rows and j=columns:

             (i-1,j-1) (i-1, j ) (i-1,j+1)

             ( i ,j-1) ( i , j ) ( i ,j+1)

             (i+1,j-1) (i+1, j ) (i+1,j+1)

    """

    # Get identities and positions of each peripheral molecule:
    # by using a 3x3 mask:

    x, y = pos
    periphery = []

    for i in range(x-1, x+2):
        for j in range(y-1, y+2):
            periphery.append(state[i][j])

    periphery = np.array(periphery).reshape(3, 3)

    """
    QUARENTINE: Method may be way too complicated to implement

    # Calculate the rate of each event occuring based
    # on a huge if else tree:

    # Set up dummy reaction constants for phenol
    # and coalin different positions:

    k_para_phenol = 0.034
    k_ortho_phenol = 0.063

    k_para_coal = 0.24
    k_ortho_coal = 0.45

    # Define the probabilities of each reaction
    # occuring (fixed for now):

    p_para_phenol = 0.63
    p_ortho_phenol = 0.37

    p_para_coal = 0.15
    p_ortho_coal = 0.85

    # Define the total chance that a crosslinking reaction
    # occurs depending on the types of peripheral molecules
    # or the chance that no reaction occurs at all (fixed for now):

    total_curing_rate_phenol = 0
    total_curing_rate_coal = 0
    total_rate_no_rxn = 0

    # If mol_type is a phenol:

    if mol_type == 1:

        # For each peripheral molecule:

        for i in range(3):

            for j in range(3):

                # if (i, j) is the central molecule, skip it:

                if i == x and j == y:

                    continue

                elif periphery[i, j] == 1:

                    # Randomly choose between k_ortho_phenol and k_para_phenol
                    # based on their respective probabilities of 
                    # occuring:

                    k = rand.choices([k_ortho_phenol, k_para_phenol], [p_ortho_phenol, p_para_phenol])

                    # Add the partial rate contribution for each
                    # peripheral molecule if is a phenol:

                    total_curing_rate_phenol += (1/8)*k

                elif periphery[i, j] == 2:
                        
                    k = rand.choices([k_ortho_coal, k_para_coal], [p_ortho_coal, p_para_coal])
                    total_curing_rate_coal += (1/8)*k # NEED EDIT

                else:

                    # If the peripheral molecule is empty 
                    # or already crosslinked, add to the 
                    # rate that no reaction occurs:

                    total_rate_no_rxn += (1/8)
                    
                    """

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

# Let's define a simpler function to get the rates
# of each event that only tracks a molecule's movement
# through the 2D space:

def get_rates_simplified(state, pos):

    """
    A function that returns the rates of the molecular
    motion of a molecule through the system in each
    of eight directions.
    
    state   : numpy array - the current state of the system in matrix
            : format.

    pos     : tuple - the position of the molecule in the matrix; this
            : is used to find the identification of the peripheral molecules
            : in the immediate vicinity of the current molecule.

            molecule i,j and its immediate peripheral molecules can be
            thought of as follows where i=rows and j=columns:

            (i-1,j-1) (i-1, j ) (i-1,j+1)

            ( i ,j-1) ( i , j ) ( i ,j+1)

            (i+1,j-1) (i+1, j ) (i+1,j+1)

    """

    # Get identities and positions of each peripheral molecule:
    # by using a 3x3 mask:

    x, y = pos
    periphery = []

    for i in range(x-1, x+2):
        for j in range(y-1, y+2):
            periphery.append(state[i][j])

    periphery = np.array(periphery).reshape(3, 3)

    # Assume the rate of the molecule moving in each
    # direction is the same for now:

    r_N, r_NE, r_E, r_SE, r_S, r_SW, r_W, r_NW = 1, 1, 1, 1, 1, 1, 1, 1
    dir_rates_matrix = [[r_NW, r_N, r_NE], 
                        [r_W,   0  , r_E], 
                        [r_SW, r_S, r_SE]]
    
    # For each peripheral molecule that is not an
    # empty space, set the rate of the molecule moving
    # to zero (central molecule cannot move):

    for i in range(3):
        for j in range(3):
            if i == 1 and j == 1:
                continue
            elif periphery[i, j] != 0:
                dir_rates_matrix[i][j] = 0

    return dir_rates_matrix

# Here is adam's code for the rate calculation
# and a wrapper for the actual KMC simulation:

def get_rates_phenol(state, pos, T):

    # Define the mol type:

    # mol_type = state[pos]

    # Get peripheral position indices and
    # identities:

    x, y = pos

    periph_pos = [(x-1, y-1), (x-1, y), (x-1, y+1),  #  0  1  2  # This object is a list
                  (x, y-1),             (x, y+1),    #  3     4  # but can be thought of
                  (x+1, y-1), (x+1, y), (x+1, y+1)]  #  5  6  7  # as a 3x3 matrix.

    periph_ident = [state[periph_pos[0]], state[periph_pos[1]], state[periph_pos[2]],
                    state[periph_pos[3]],                       state[periph_pos[4]],
                    state[periph_pos[5]], state[periph_pos[6]], state[periph_pos[7]]]

    # Assign easy names to peripheral position identities:
    """
    top_left, top, top_right = periph_ident[0], periph_ident[1], periph_ident[2]
    left, right = periph_ident[3], periph_ident[4]
    bot_left, bot, bot_right = periph_ident[5], periph_ident[6], periph_ident[7]
    """
    # Tally the identities of each peripheral molecule
    # (this will be updated to accommodate ortho or para
    # positions and pre-cured molecules):

    n_phenol = periph_ident.count(1)
    n_coal = periph_ident.count(2)
    # n_void = periph_ident.count(0)

    # Define rxn rate constants
    # (simplified for now):

    kpp = 0.8 * T
    kpc = 0.2 * T

    # Define the KMC rates of each event (weighted by
    # the number of peripheral molecules of each type):

    r_kmc_pp = (kpp - (kpc * 0.2)) * (n_phenol)
    r_kmc_pc = (kpc - (kpp * 0.2)) * (n_coal)
    r_kmc_no_rxn = (0 if (r_kmc_pp or r_kmc_pc) else 1)

    if r_kmc_no_rxn < 0:
        r_kmc_no_rxn = 0

    # Define the KMC rates of each movement event:

    if r_kmc_pp or r_kmc_pc:
        move_prob = 1 / (n_phenol + n_coal)
    else:
        move_prob = 1

    r_kmc_move_ul = (1 if periph_ident[0] == 0 else 0) * move_prob
    r_kmc_move_u = (1 if periph_ident[1] == 0 else 0) * move_prob
    r_kmc_move_ur = (1 if periph_ident[2] == 0 else 0) * move_prob
    r_kmc_move_l = (1 if periph_ident[3] == 0 else 0) * move_prob
    r_kmc_move_r = (1 if periph_ident[4] == 0 else 0) * move_prob
    r_kmc_move_dl = (1 if periph_ident[5] == 0 else 0) * move_prob
    r_kmc_move_d = (1 if periph_ident[6] == 0 else 0) * move_prob
    r_kmc_move_dr = (1 if periph_ident[7] == 0 else 0) * move_prob

    return (r_kmc_pp, r_kmc_pc, r_kmc_no_rxn,
            r_kmc_move_ul, r_kmc_move_u, r_kmc_move_ur,
            r_kmc_move_l,                r_kmc_move_r,
            r_kmc_move_dl, r_kmc_move_d, r_kmc_move_dr,
            periph_pos, periph_ident)                       # Peripheral positions and identities

def get_rates_coal(state, pos, T):

    # Define the mol type:

    # mol_type = state[pos]

    # Get peripheral position indices and
    # identities:

    x, y = pos

    periph_pos = [(x-1, y-1), (x-1, y), (x-1, y+1),  #  0  1  2  # This object is a list
                  (x, y-1),             (x, y+1),    #  3     4  # but can be thought of
                  (x+1, y-1), (x+1, y), (x+1, y+1)]  #  5  6  7  # as a 3x3 matrix.

    periph_ident = [state[periph_pos[0]], state[periph_pos[1]], state[periph_pos[2]],
                    state[periph_pos[3]],                       state[periph_pos[4]],
                    state[periph_pos[5]], state[periph_pos[6]], state[periph_pos[7]]]

    # Assign easy names to peripheral position identities:
    """
    top_left, top, top_right = periph_ident[0], periph_ident[1], periph_ident[2]
    left, right = periph_ident[3], periph_ident[4]
    bot_left, bot, bot_right = periph_ident[5], periph_ident[6], periph_ident[7]
    """
    # Tally the identities of each peripheral molecule
    # (this will be updated to accommodate ortho or para
    # positions and pre-cured molecules):

    n_phenol = periph_ident.count(1)
    n_coal = periph_ident.count(2)
    # n_void = periph_ident.count(0)

    # Define rxn rate constants
    # (simplified for now)
    # coal cannot react with itself, so we only need
    # to define the rate of coal reacting with phenol:

    kcp = 0.2 * T # Again we need to change this

    # Define the KMC rates of each event (weighted by
    # the number of peripheral molecules of each type):

    r_kmc_cp = kcp * (n_coal)
    r_kmc_no_rxn = (0 if kcp else 1)

    if r_kmc_no_rxn < 0:
        r_kmc_no_rxn = 0

    # Define the KMC rates of each movement event:

    if r_kmc_cp:
        move_prob = 1 / (n_phenol + n_coal)
    else:
        move_prob = 1

    r_kmc_move_ul = (1 if periph_ident[0] == 0 else 0) * move_prob
    r_kmc_move_u = (1 if periph_ident[1] == 0 else 0) * move_prob
    r_kmc_move_ur = (1 if periph_ident[2] == 0 else 0) * move_prob
    r_kmc_move_l = (1 if periph_ident[3] == 0 else 0) * move_prob
    r_kmc_move_r = (1 if periph_ident[4] == 0 else 0) * move_prob
    r_kmc_move_dl = (1 if periph_ident[5] == 0 else 0) * move_prob
    r_kmc_move_d = (1 if periph_ident[6] == 0 else 0) * move_prob
    r_kmc_move_dr = (1 if periph_ident[7] == 0 else 0) * move_prob

    return (r_kmc_cp, r_kmc_no_rxn,
            r_kmc_move_ul, r_kmc_move_u, r_kmc_move_ur,
            r_kmc_move_l,                r_kmc_move_r,
            r_kmc_move_dl, r_kmc_move_d, r_kmc_move_dr,
            periph_pos, periph_ident)

# Now we need a function to pick the event that occurs 
# based on the precalculated rates:

def choose_event(rates):

    """
    Chooses the event to occur based on the rates in
    tuple format.
    
    """

    # Sum the rates:

    total_rate = sum(rates)
    cumulative_rate = np.cumsum(rates)

    # Randomly choose an event based on the rates:

    choice = rand.uniform(0, total_rate)

    if 0 <= choice < cumulative_rate[0]:
        return "pp_rxn"
    elif cumulative_rate[0] <= choice < cumulative_rate[1]:
        return "pc_rxn"
    elif cumulative_rate[1] <= choice < cumulative_rate[2]:
        return "no_rxn"
    elif cumulative_rate[2] <= choice < cumulative_rate[3]:
        return "move_ul"
    elif cumulative_rate[3] <= choice < cumulative_rate[4]:
        return "move_u"
    elif cumulative_rate[4] <= choice < cumulative_rate[5]:
        return "move_ur"
    elif cumulative_rate[5] <= choice < cumulative_rate[6]:
        return "move_l"
    elif cumulative_rate[6] <= choice < cumulative_rate[7]:
        return "move_r"
    elif cumulative_rate[7] <= choice < cumulative_rate[8]:
        return "move_dl"
    elif cumulative_rate[8] <= choice < cumulative_rate[9]:
        return "move_d"
    elif cumulative_rate[9] <= choice < cumulative_rate[10]:
        return "move_dr"

# Adam's rate calculator still needs to be updated
# but is our best bet at the moment. Here is a fucntion
# that calculates the new state of the system:

def get_new_state(current_state, T):
    
    """
    A function that calculates the new state of the system
    after a KMC event occurs.

    current_state   : numpy array - the current state of the system in matrix
                    : format.

    T               : float - the temperature of the system.

    """
    
    # Get the dimensions of the current state:

    sim_size = current_state.shape[0]

    # Find non-void positions in the current state
    # that are not around the edges of the system:

    phenol_pos = np.where(current_state == 1)
    coal_pos = np.where(current_state == 2)

    phenol_pos = list(zip(phenol_pos[0], phenol_pos[1]))
    coal_pos = list(zip(coal_pos[0], coal_pos[1]))

    # Only consider the non-void positions that are
    # not around the edges of the system:

    phenol_pos = [pos for pos in phenol_pos if pos[0] != 0 and pos[0] != sim_size-1 and pos[1] != 0 and pos[1] != sim_size-1]
    coal_pos = [pos for pos in coal_pos if pos[0] != 0 and pos[0] != sim_size-1 and pos[1] != 0 and pos[1] != sim_size-1]

    # For each non-void position, calculate the rates
    # and choose event that occurs: