"""This module aligns and colours dataframes as compared to noptm using 7stz as reference """

import numpy as np
from numpy import array,dot
from MDAnalysis.analysis import align
import pandas as pd
from functions.CIF_Functions import make_df_atomlist,hetatm_df,make_df_new
from functions.Sliding_Functions import local_rmsd_plotter

def align_by_subsection(df_a,df_b,rel_sub_start:int,rel_sub_end:int):
    """Function aligns coordinate arrays or dataframes over a subset  
    
    Args:
        Equally indexed dataframes with coords in columns [10:13]
        or 3d array 

    Returns:
        aligned protein A, B, aligned region RMSD, Rotational matrix, subset-Mean-A, Subset-Mean-B 
    """
    # type changing so either dataframes or coord_arrays are allowed
    if isinstance(df_a,pd.DataFrame):
        x = array(df_a[df_a.columns[10:13]]).astype(float)
    else:
        x = df_a

    if isinstance(df_b,pd.DataFrame):
        y = array(df_b[df_b.columns[10:13]]).astype(float)
    else:
        y= df_b

     ## define subsections
    x_s = x[rel_sub_start:rel_sub_end]
    y_s = y[rel_sub_start:rel_sub_end]

    # Mean centering full protein by subset
    x_f = x - (sum(x_s)/len(x_s))
    y_f = y - (sum(y_s)/len(y_s))

    mean_centre_a = (sum(x_s)/len(x_s))
    mean_centre_b = (sum(y_s)/len(y_s))

    ## Aligning these centered subsets
    rot_matrix, rms = align.rotation_matrix(x_f[rel_sub_start:rel_sub_end],
                                            y_f[rel_sub_start:rel_sub_end])
    # instead of rotation B we rotate A with the transpose of R
    x_f = dot(x_f,np.transpose(rot_matrix))

    return array(x_f),array(y_f),rms,rot_matrix,mean_centre_a,mean_centre_b


def align_color_cif(CIF_PATH,CHAIN):
    """Takes crystal, and noptm crystal align them to a reference (7stz) and colour  
    
    Args:
        Dataframe with coords in columns [10:13] to be aligned to noptm

    Returns:
        aligned to noptm and coloured by ref 7stz vs noptm
    """
# Read dataframes
    atomlist = ["CA","N","C"]
    ptm_df = make_df_atomlist(CIF_PATH,atomlist)
    noptm_df = make_df_atomlist("cif_files/Predictions/SEEDMATCHED/fold_adnan_seed1_noptms/fold_adnan_seed1_noptms_model_1.cif",atomlist)
    crystal_df = make_df_atomlist("cif_files/Crystals/7stz (1).cif",atomlist)

# Find out how to align (this gives R, a rotational matrix TO BE USED ON B)
    # Select the first 2 domains (roughly) (speed increase for alignment)
    STOP = 600 # arbitrary
    df_a,df_b = ptm_df[:STOP],noptm_df[:STOP]
    align_start = 72*3
    align_end   = 115*3  # needs to be relative to the start of the full protein
    _,_,_,R,mcA,mcB = align_by_subsection(df_a,df_b,align_start,align_end)
    # this gives us the rotational and tranlational information to apply to the full molecule later in the function


# Coloring

    # Do sliding window to get significant difference peeks
    win = 25*len(atomlist) # 25 residues
    _,y1,_ = local_rmsd_plotter(ptm_df,crystal_df,stepsize=1,win_size=win)
    _,y2,_ = local_rmsd_plotter(noptm_df,crystal_df,stepsize=1,win_size=win)
    # Plot the difference
    data = array(y2)-array(y1)
    # Remove outliers formed by the 5 initial values at the beginning of prediction and make the colourlist
    data_truncated = [i for i in data if i>-.5 and i<.5]
    sym_neg_data = [i for i in data_truncated if i<0] + [-1*i for i in data_truncated if i<0]
    sig3 = 3*np.std(sym_neg_data)
    sig2 = 2*np.std(sym_neg_data)
    colorlist= ["red" if i> sig3 or i< -sig3 else "orange" if i> sig2 or i< -sig2 else "blue" for i in data]
    # The sliding window function starts its value halfway in the window so the first and last 37 atoms are informationsless, so we color them insignificant
    colorlist = ["blue" for _ in range(37)] + colorlist + ["blue" for _ in range(37)]

# Getting ALL atoms (also sidechains)
    #full_ptm_df,sequenceA = make_df_new("cif_files/Predictions/SEEDMATCHED/fold_adnan_seed42_ptms/fold_adnan_seed42_ptms_model_3.cif")
    full_ptm_df,_ = make_df_new(CIF_PATH)
    #full_noptm_df,_ = make_df_new("cif_files/Predictions/SEEDMATCHED/fold_adnan_seed1_noptms/fold_adnan_seed1_noptms_model_1.cif")
 
    # Read xyz
    f_ptm_df = full_ptm_df[0] 
    #f_noptm_df = full_noptm_df[0]
    tb_aligned_a = array(f_ptm_df[f_ptm_df.columns[10:13]]).astype(float)
    #tb_aligned_b = array(f_noptm_df[f_noptm_df.columns[10:13]]).astype(float)

    # mean center on the aligned part
    mc_a = tb_aligned_a - mcA
    #mc_b = tb_aligned_b -  mcB
    # run the alignment 
    aligned_a = dot(mc_a,np.transpose(R))

    # the same for HETATMS
    het_df = hetatm_df(CIF_PATH)
    het_xyz = array(het_df[het_df.columns[10:13]]).astype(float) # select xyz values
    het_xyz = het_xyz - mcA # mean-center
    het_xyz = dot(het_xyz,np.transpose(R)) # do the rotation


    # setting axes to aligned values
    het_df[["Cartn_x","Cartn_y","Cartn_z"]] = het_xyz
    het_df[["label_asym_id","auth_asym_id"]] = CHAIN
    het_df["B_iso_or_equiv"] = 0
# Column B factor changing

    colorlist_c = [3 if i=="red" else 2 if i=="orange" else 1 for i in colorlist]

    #f_ptm_df.loc[f_ptm_df["label_atom_id"]=="CA","B_iso_or_equiv"] = [colorlist_c[i] for i in range(len(colorlist)) if i%3==1]
    #f_ptm_df.loc[f_ptm_df["label_atom_id"].isin(["CA","C","N"]), "B_iso_or_equiv"] = colorlist_c
    #f_ptm_df.loc[f_ptm_df["label_atom_id"].isin(["CA","C","N"]), "occupancy"] = colorlist_c
    CA_only_SIGN = [colorlist_c[i] for i in range(len(colorlist)) if i%3==1]
    for seq_id in range(1,541):
        f_ptm_df.loc[f_ptm_df["label_seq_id"]==f'{seq_id}',"B_iso_or_equiv"] = CA_only_SIGN[seq_id-1]

    # The other atoms

    #f_ptm_df.loc[~f_ptm_df["label_atom_id"].isin(["CA","C","N"]), "B_iso_or_equiv"] = 0

# changing xyz's to aligned xyz
    # Modify the first dataframe to hold the transformed values
    copy  = f_ptm_df
    # setting axes to aligned values
    copy[["Cartn_x","Cartn_y","Cartn_z"]] = aligned_a
    copy[["label_asym_id","auth_asym_id"]] = CHAIN
    final_df = pd.concat([copy,het_df])
    return final_df

# Crystal Files only need to be aligned to noptms and not colored to a ref

def align_to_noptms(CIF_PATH,CHAIN):
    """Takes crystal, and aligns to noptm crystal  
    
    Args:
        Dataframe with coords in columns [10:13] to be aligned to noptm

    Returns:
        aligned to noptm 
    """
# Read dataframes
    atomlist = ["CA","N","C"]
    ptm_df = make_df_atomlist(CIF_PATH,atomlist)
    noptm_df = make_df_atomlist("cif_files/Predictions/SEEDMATCHED/fold_adnan_seed1_noptms/fold_adnan_seed1_noptms_model_1.cif",atomlist)

# Find out how to align (this gives R, a rotational matrix TO BE USED ON B)
    # Select the first 2 domains (roughly) (speed increase for alignment)
    STOP = 600 # arbitrary
    df_a,df_b = ptm_df[:STOP],noptm_df[:STOP]
    align_start = 72*3
    align_end   = 115*3  # needs to be relative to the start of the full protein
    _,_,_,R,mcA,mcB = align_by_subsection(df_a,df_b,align_start,align_end) 
    # this gives us the rotational and tranlational information to apply to the full molecule later in the function


# Getting ALL atoms (also sidechains)
    #full_ptm_df,sequenceA = make_df_new("cif_files/Predictions/SEEDMATCHED/fold_adnan_seed42_ptms/fold_adnan_seed42_ptms_model_3.cif")
    full_ptm_df,_ = make_df_new(CIF_PATH)
    # Read xyz
    f_ptm_df = full_ptm_df[0] 
    tb_aligned_a = array(f_ptm_df[f_ptm_df.columns[10:13]]).astype(float)
    # mean center on the aligned part
    mc_a = tb_aligned_a - mcA
    # run the alignment 
    aligned_a = dot(mc_a,np.transpose(R))

# Column B factor changing
    atomlist=["CA","C","N"]
    f_ptm_df.loc[f_ptm_df["label_atom_id"].isin(atomlist), "B_iso_or_equiv"] = 1
    # The other atoms
    f_ptm_df.loc[~f_ptm_df["label_atom_id"].isin(atomlist), "B_iso_or_equiv"] = 0

# changing xyz's to aligned xyz
    # Modify the first dataframe to hold the transformed values
    copy  = f_ptm_df
    # setting axes to aligned values
    copy[["Cartn_x","Cartn_y","Cartn_z"]] = aligned_a
    copy[["label_asym_id","auth_asym_id"]] = CHAIN
    return copy
