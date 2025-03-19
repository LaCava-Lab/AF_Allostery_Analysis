import numpy as np
import pandas as pd
from numpy import array,dot
from MDAnalysis.analysis import align  
from MDAnalysis.analysis.rms import rmsd

def local_rmsd_plotter(df_a :pd.DataFrame,
                       df_b :pd.DataFrame,
                       stepsize=2,
                       win_size=10) -> tuple[list[int],list[int],list[int]]:
    """Function calculate aligned RMSD over a set window 
    Args:
        equally indexed dataframes with coords in columns [10:13]

    Returns:
        list of residue indexes , list of RMSD's
    """
    ## for stepsizeChecking
    #if stepsize not in find_step_size(len(df_a),win_size):
    #    print(f"Inconsistent length of window, pick {find_step_size(len(df_a),win_size)}")
    #    return 0,0
    assert len(df_a) == len(df_b), "Sequences are differing length"
     
    win_start = 0
    win_end = 0
    score_list = []
    rmsd_list = []
    start_list=[]

    #for _ in range(round(len(df_b)/stepsize)):
    while win_end < len(df_a):
        # Keep track of the residue in the middle of the window
        start_list.append(win_start+int(win_size/2))
        # Define window
        win_end = win_start + win_size

        # Do inside that window

        # Fetch average pLDDT
        score_list.append(np.mean(df_a["B_iso_or_equiv"][win_start:win_end].astype(float))) 
        # get a rmsd comparing that specific region
        x= array(df_a[df_a.columns[10:13]][win_start:win_end]).astype(float)
        x_mean_centered = x - (sum(x)/len(x))

        y= array(df_b[df_b.columns[10:13]][win_start:win_end]).astype(float)
        y_mean_centered = y - (sum(y)/len(y))

        _, rms = align.rotation_matrix(x_mean_centered, y_mean_centered)
        rmsd_list.append(rms)

        # Define next window
        win_start+=stepsize
        
    return start_list,rmsd_list,score_list


def global_rmsd_plotter(df_a,
                        df_b,
                        stepsize=2,
                        win_size=10) -> tuple[list[pd.DataFrame],list[pd.DataFrame],float]:
    """Function calculate RMSD over a set window aligned globally
    Args:
        Equally indexed dataframes with coords in columns [10:13]

    Returns:
        list of residue indexes , list of RMSD's
    """
     # type changing
    if type(df_a)==pd.DataFrame:
        x = array(df_a[df_a.columns[10:13]]).astype(float)
    else:
        x = df_a
    
    if type(df_b) ==pd.DataFrame:
        y = array(df_b[df_b.columns[10:13]]).astype(float)
    else:
        y= df_b


    # First Global allignment
    # Define mean centered Prediction coordinates
    #x = array(df_a[df_a.columns[10:13]]).astype(float)
    x_mean_centered = x - (sum(x)/len(x))

    # Define mean centered Crystal coordinates
    #y = array(df_b[df_b.columns[10:13]]).astype(float)
    y_mean_centered = y - (sum(y)/len(y))

    # align the 2, gives global rmsd, and rotational matrix
    R, rms = align.rotation_matrix(x_mean_centered, y_mean_centered)
    y_mean_centered_alligned = dot(y_mean_centered,R)
    grms = rms
    #########################################################
    win_end = 0
    win_start = 0
    gbl_rmsd_list=[]
    start_list=[]

    while win_end< len(df_a):
    #for i in range(round((len(df_b))/stepsize)):
        start_list.append(win_start)
        # Define window
        win_end = win_start + win_size
        af_coords = x_mean_centered[win_start:win_end]
        xr_coords = y_mean_centered_alligned[win_start:win_end]
        # get the RMSD for that window WITHOUT ALLIGNMENT
        rms = rmsd(af_coords,xr_coords)
        gbl_rmsd_list.append(rms)
        # Define next window
        win_start+=stepsize
    return start_list,gbl_rmsd_list,grms

from MDAnalysis.analysis.rms import rmsd

def local_rmsd_plotter_no_aligning(array_a,
                       array_b,
                       stepsize=2,
                       win_size=10) -> tuple[list[int],list[int]]:
    """Function calculate aligned RMSD over a set window 
    Args:
        equally indexed dataframes with coords in columns [10:13]

    Returns:
        list of residue indexes , list of RMSD's
    """
    ## for stepsizeChecking
    #if stepsize not in find_step_size(len(df_a),win_size):
    #    print(f"Inconsistent length of window, pick {find_step_size(len(df_a),win_size)}")
    #    return 0,0
    assert len(array_a) == len(array_b), "Sequences are differing length"

    win_start = 0
    win_end = 0
    rmsd_list = []
    start_list=[]

    #for _ in range(round(len(df_b)/stepsize)):
    while win_end < len(array_a):
        # Keep track of the residue in the middle of the window
        start_list.append(win_start+int(win_size/2))
        # Define window
        win_end = win_start + win_size

        # Do inside that window

        # get a rmsd comparing that specific region
        x= array(array_a[win_start:win_end]).astype(float)
        y= array(array_b[win_start:win_end]).astype(float)

        rmsd_list.append(rmsd(x,y))

        # Define next window
        win_start+=stepsize
 
    return start_list,rmsd_list

def get_significance(x,y1,y2):
    """Function calculates sigmas and colours accordingly
    Args:
       rmsd lists y1 and y2 to take a delta from 

    Returns:
        colourlist of length data
    """
    data = array(y2)-array(y1)
    # Remove outliers formed by the 5 initial values at the beginning of prediction
    data_truncated = [i for i in data if i>-.5 and i<.5]
    sym_neg_data = [i for i in data_truncated if i<0] + [-1*i for i in data_truncated if i<0]

    # Plot this symetrical density plot
    sig3 = 3*np.std(sym_neg_data)
    sig2 = 2*np.std(sym_neg_data)

    return ["red" if i> sig3 or i< -sig3 else "orange" if i> sig2 or i< -sig2 else "b" for i in data]
