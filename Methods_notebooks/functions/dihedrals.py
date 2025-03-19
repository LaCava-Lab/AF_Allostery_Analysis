from math import sqrt,atan2,degrees

def dot_product(v1, v2):
    """ Calculate the dot product of two vectors """
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
#print(dot_product([1, 2, 3], [1, 3, 2]))

def cross_product(v1, v2):
    """ Calculate the cross product of two vectors """
    i = v1[1]*v2[2] - v1[2]*v2[1]
    j = v1[2]*v2[0] - v1[0]*v2[2]
    k = v1[0]*v2[1] - v1[1]*v2[0]
    return [i,j,k]
#print(cross_product([1, 2, 3], [1, 3, 2]))

def magnitude(v):
    """ Calculate the size of a vector """
    return sqrt(v[0]**2 + v[1]**2 + v[2]**2)


def read_df(dataframe):
    # dictionaries to store the output
    # pdb atom coordinates:
    #     pdbcoord[chain][residue_number][atom_type] = coordinates
    pdbcoord = {}
    # residue type per chain and residue number (i.e. store the sequence)
    #     pdbseq[chain][resnum] = restype
    pdbseq = {}
    columnnumber = list(dataframe.columns).index("auth_seq_id")
    for row in dataframe.values:
        atom_type = row[3]
        # AMINO ACID type (e.g. alanine)
        aa_type = row[5]
        # residue number

        res_num = int(row[columnnumber])
        # Protein chain
        chain = row[6]
        # coordinates
        xcoord = float(row[10])
        ycoord = float(row[11])
        zcoord = float(row[12])

        if not chain in pdbcoord:
            pdbcoord[chain] = {}
            pdbseq[chain] = {}
        # if resnum does not exists create new entry
        if not res_num in pdbcoord[chain]:
            pdbcoord[chain][res_num] = {}

        # store coordinates as a vector
        pdbcoord[chain][res_num][atom_type] = [xcoord,ycoord,zcoord]
        # store sequence
        pdbseq[chain][res_num] = aa_type
    return pdbcoord,pdbseq


def calculateDihedral(a1, a2, a3, a4):
    """ Calculates the normal vector of the planes
    defined by four atom coordinates """

    ### START CODING HERE
    # calculate normal vectors to the planes defined by a1,a2,a3 and a2,a3,a4
    # you may use the functions "cross_product","dot_product" and "magnitude" defined above
    # you can also use the python math function "atan2" and "degrees"

    # Calculate connecting vectors (bonds)
    b1 = [a2[i] - a1[i] for i in range(len(a1))]
    b2 = [a3[i] - a2[i] for i in range(len(a2))]
    b3 = [a4[i] - a3[i] for i in range(len(a3))]

    # Calculate normal vectors to the planes
    n_p1 = cross_product(b1, b2)
    n_p2 = cross_product(b2, b3)

    #calculate a sign
    n_p1_normal = [a/magnitude(n_p1) for a in n_p1 ]
    n_p2_normal = [a/magnitude(n_p2) for a in n_p2 ]

    ortho_vector = cross_product(n_p1_normal,n_p2_normal)

    sign = dot_product(ortho_vector,b2) / magnitude(b2)

    if sign >= 0:
        sign = 1
    else:
        sign = -1

    # Calculate the dihedral angle
    sin = magnitude(ortho_vector) * sign
    cos = dot_product(n_p1_normal, n_p2_normal)
    
    dihedral = degrees(atan2(sin, cos))
    ### END CODING HERE
    return dihedral

def assign_dihedrals(pdbcoord):
    phi_list = []
    psi_list = []
    for chain in pdbcoord.keys():
        list_residue_numbers = sorted(pdbcoord[chain].keys())
        for res_num in sorted(pdbcoord[chain].keys()):
            # if certain residues are missing in the PDB file, you will
            # get a KeyError. Make sure your program does not crash, but
            # gives a warning when this happens
            try:
                phi = None
                psi = None
            
                if res_num > 1:
                    phi = calculateDihedral(pdbcoord[chain][res_num-1]['C'],
                                        pdbcoord[chain][res_num]['N'],
                                        pdbcoord[chain][res_num]['CA'],
                                        pdbcoord[chain][res_num]['C'])
                    phi_list.append(phi)
                    
                if res_num< max(list_residue_numbers):
                    psi = calculateDihedral(pdbcoord[chain][res_num]['N'],
                                        pdbcoord[chain][res_num]['CA'],
                                        pdbcoord[chain][res_num]['C'],
                                        pdbcoord[chain][res_num+1]['N'])
                    psi_list.append(psi)
            except KeyError:
                    print('WARNING: KeyError:', KeyError, 'in residue', chain, res_num)  
    return phi_list,psi_list