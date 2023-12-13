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

    kpp = 0.8 * T # Use Jessica's rate functions!
    kpc = 0.2 * T

    # Define the KMC rates of each event (weighted by
    # the number of peripheral molecules of each type):

    r_kmc_pp = (kpp - (kpc * 0.2)) * (n_phenol) # To be changed!
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

def choose_event(rates, mol_type):

    """
    Chooses the event to occur based on the rates in
    tuple format.
    
    """

    # Sum the rates:

    total_rate = sum(rates)
    cumulative_rate = np.cumsum(rates)

    # Randomly choose an event based on the rates:

    choice = rand.uniform(0, total_rate)
    
    if mol_type == 1:

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
        
    elif mol_type == 2:

        if 0 <= choice < cumulative_rate[0]:
            return "pc_rxn"
        elif cumulative_rate[0] <= choice < cumulative_rate[1]:
            return "no_rxn"
        elif cumulative_rate[1] <= choice < cumulative_rate[2]:
            return "move_ul"
        elif cumulative_rate[2] <= choice < cumulative_rate[3]:
            return "move_u"
        elif cumulative_rate[3] <= choice < cumulative_rate[4]:
            return "move_ur"
        elif cumulative_rate[4] <= choice < cumulative_rate[5]:
            return "move_l"
        elif cumulative_rate[5] <= choice < cumulative_rate[6]:
            return "move_r"
        elif cumulative_rate[6] <= choice < cumulative_rate[7]:
            return "move_dl"
        elif cumulative_rate[7] <= choice < cumulative_rate[8]:
            return "move_d"
        elif cumulative_rate[8] <= choice < cumulative_rate[9]:
            return "move_dr"
        
    else:

        raise("Invalid molecule type in choose_event function.")

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

    # Create a copy of the current state to be used
    # as a blank canvas for the new state:

    new_state = np.copy(current_state)

    # Find non-void positions in the current state
    # that are not around the edges of the system:

    phenol_pos = np.where(current_state == 1)
    coal_pos = np.where(current_state == 2)

    phenol_pos = list(zip(phenol_pos[0], phenol_pos[1]))
    coal_pos = list(zip(coal_pos[0], coal_pos[1]))

    # Only consider the molecular positions that are
    # not around the edges of the system:

    phenol_pos = [pos for pos in phenol_pos if pos[0] != 0 and pos[0] != sim_size-1 and pos[1] != 0 and pos[1] != sim_size-1]
    coal_pos = [pos for pos in coal_pos if pos[0] != 0 and pos[0] != sim_size-1 and pos[1] != 0 and pos[1] != sim_size-1]

    # For each molecular position, calculate the rates
    # and choose event that occurs:

    all_mols = phenol_pos + coal_pos

    # What are we extracting after the new state is calculated?
    # We need to extract the final heat of reaction from summing
    # all of the reactions that occur in the system (phenol to
    # phenol, and phenol to coal), how many coal particles reacted
    # with phenol, and the degree of crosslinking
    # in the system (number of crosslinked phenol molecules). We
    # will need to set these variables up here:

    # Number of crosslinks occuring in this iteration of the system:

    n_crosslinks = 0

    # Number of coal particles that reacted with phenol in this
    # iteration of the system:

    n_coal_rxn = 0

    # Heat of reaction of this iteration of the system:

    mol_heat_rxn = 0 # Jessica has a function for this!

    for coords in all_mols:

        # Technically, this would overwrite some of the changes made
        # previously to the new state if a molecule is in the same
        # 3x3 matrix as another molecule since the frame of reference
        # can overlap, but we may not have time to change this.

        # Get the rates of each event occuring (since the get_rates
        # functions output the rates and the peripheral positions, 
        # we need to save both data separately):

        if current_state[coords] == 1:

            all_rate_data = get_rates_phenol(current_state, coords, T)
            rates = all_rate_data[0:11]
            periph_pos, periph_ident = all_rate_data[11], all_rate_data[12]
            
            # Choose the event that occurs:

            event = choose_event(rates, 1)

        elif current_state[coords] == 2:

            all_rate_data = get_rates_coal(current_state, coords, T)
            rates = all_rate_data[0:10]
            periph_pos, periph_ident = all_rate_data[10], all_rate_data[11]

            # Choose the event that occurs:

            event = choose_event(rates, 2)

        else:

            # If the current state is not a phenol or coal (i.e. it
            # has crosslinked or reacted with coal), then we do not
            # need to calculate the rates of the events occuring:

            continue

        # Update the new state based on the event that
        # occurs:

        if event == "pp_rxn":

            # Phenol molecule becomes a 3 to represent that 
            # it has crosslinked:

            new_state[coords] = 3

            # However, the molecule which it reacted with also
            # has to become a 3 since both molecules crosslinked.
            # Additionally, since any of the peripheral phenol
            # molecules (if multiple exist) could have reacted
            # with an equal probability, we need to randomly
            # determine which peripheral phenol molecule reacted
            # and update that new molecule to a 3 as well:

            periph_phenol_index = np.where(np.array(periph_ident) == 1)[0]

            # Randomly choose a peripheral phenol molecule:

            random_phenol_choice = rand.choice(periph_phenol_index)

            # Update the new state:

            new_state[periph_pos[random_phenol_choice]] = 3

            # Update the number of crosslinks and the heat of
            # reaction:

            n_crosslinks += 1
            mol_heat_rxn += 0.00043 # This is a placeholder value

        elif event == "pc_rxn":

            # Phenol molecule becomes a 4 to represent that
            # it has reacted with coal:

            new_state[coords] = 4

            # Find peripheral coal molecules:

            periph_coal_index = np.where(np.array(periph_ident) == 2)[0]

            # Randomly choose a peripheral coal molecule:

            random_coal_choice = rand.choice(periph_coal_index)

            # Update the new state:

            new_state[periph_pos[random_coal_choice]] = 4

            # Update the number of coal particles that reacted:

            n_coal_rxn += 1
            mol_heat_rxn += 0.00015 # This is a placeholder value

        elif event == "no_rxn":
            pass
        elif event == "move_ul":
            new_state[coords] = 0
            new_state[periph_pos[0]] = current_state[coords]
        elif event == "move_u":
            new_state[coords] = 0
            new_state[periph_pos[1]] = current_state[coords]
        elif event == "move_ur":
            new_state[coords] = 0
            new_state[periph_pos[2]] = current_state[coords]
        elif event == "move_l":
            new_state[coords] = 0
            new_state[periph_pos[3]] = current_state[coords]
        elif event == "move_r":
            new_state[coords] = 0
            new_state[periph_pos[4]] = current_state[coords]
        elif event == "move_dl":
            new_state[coords] = 0
            new_state[periph_pos[5]] = current_state[coords]
        elif event == "move_d":
            new_state[coords] = 0
            new_state[periph_pos[6]] = current_state[coords]
        elif event == "move_dr":
            new_state[coords] = 0
            new_state[periph_pos[7]] = current_state[coords]
    
    # We need to update the resulting temperature of the system
    # based on the output heat of reaction:
    # Assume: delta H = m * c * (T_final - T_initial)
    # So: T_final = T_initial + (delta H / (m * c))
    
    # Assume a constant heat capacity of resin and coal:

    heat_capacity = 1 # This is a placeholder value
    m = 1 # This is a placeholder value

    # Update the temperature of the system:

    T += mol_heat_rxn / (m * heat_capacity)

    # Jessica has a reaction to accurately define this!

    # Return the new state of the system and the desired
    # output variables:

    return new_state, T, n_crosslinks, n_coal_rxn, mol_heat_rxn

# Now we need to wrap the get_new_state function in a
# function that will run the simulation for a given
# number of system state iterations:

def resin_cure_simulation(n, ratio, T, n_iter):

    """
    A function that runs the simulation for a given
    number of system state iterations.

    n       : int - the size of the system in the x
            : and y dimensions.

    ratio   : list - the ratio of voids to phenol to
            : coal in the system at time 0.

    T       : float - the initial temperature of the
            : system.

    n_iter  : int - the number of iterations of the
            : system state.

    """

    # Generate the initial state of the system:

    state = initial_random_matrix(n, ratio)

    # Initialize the output variables:

    state_list = [state]

    crosslinks = []
    coal_rxn = []
    heat_rxn = []
    temps = [T]

    # Run the simulation for the desired number of
    # iterations:

    for i in range(n_iter):

        # Calculate the new state of the system:

        new_state, T, n_crosslinks, n_coal_rxn, mol_heat_rxn = get_new_state(state, T)

        # Update the state of the system:

        state_list.append(new_state)
        state = new_state

        # Append the output variables:

        crosslinks.append(n_crosslinks)
        coal_rxn.append(n_coal_rxn)
        heat_rxn.append(mol_heat_rxn)
        temps.append(T)

    # Return the output variables:

    return state_list, temps, crosslinks, coal_rxn, heat_rxn