
import numpy as np

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