"""Functions to simplify and quicken the setup process

    select_ca()
    make_plot()?  
"""
import math
import pandas as pd
from functions.CIF_Functions import make_df_new,hetatm_df
from numpy import array,dot
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd


def dist(coord1, coord2) -> float:
    """Takes 2 xyz coordinates and calculates their distance

    Args:
        list[x,y,z]

    Returns:
        float
    """
    x1, y1, z1 = coord1
    x2, y2, z2 = coord2
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

def select_ca(raw_df)->pd.DataFrame:
    """Takes in a df of atoms and returns only the C-Alpha atoms,
    also removes any alternative residues (bugfix441)

    Args:
        pd.DataFrame of all atoms

    Returns:
        pd.DataFrame of C-alpha atoms
    """
    
    im_df = raw_df[raw_df["label_atom_id"]=="CA"]
    df_ca = im_df[(im_df['label_alt_id'] == ".")|(im_df['label_alt_id']=="A")] #We could average A/B
    return df_ca


def get_hetatm_rmsd(df_a_path, df_b_path)->float:
    """Takes in a path of cif file, makes of ATOM and HETATM
    dataframes and return the rmsd betweenn ordered Calcium ions
    
    Args:
        PATH to cif file containing ATOM and HETATM lines

    Returns:
        RMSD
    """
    
    # Making dataframes
    # The make df new function makes a list of sequences and dataframes.
    # we select the first item from the list of dataframes
    df_a = make_df_new(df_a_path)[0][0]
    df_b = make_df_new(df_b_path)[0][0]
    # The hetatom_df function gets all HETATM lines 
    # selects CA (this is the same function as C-alpha but the identifier is the same "CA")
    ac = hetatm_df(df_a_path)
    bc = hetatm_df(df_b_path)
    df_a = select_ca(df_a)
    df_b = select_ca(df_b)
    ac = select_ca(ac)
    bc = select_ca(bc)
    
    if "Crystals" in df_a_path:
        ac = ac[ac["auth_asym_id"]=="C"]
    elif "Crystals" in df_b_path:
        bc = bc[bc["auth_asym_id"]=="C"]
    
    ### Aligning
    # Mean Centering
    x= array(df_a[df_a.columns[10:13]]).astype(float)
    x_mean_centered = x - (sum(x)/len(x))

    y= array(df_b[df_b.columns[10:13]]).astype(float)
    y_mean_centered = y - (sum(y)/len(y))

    # mean center Calcium
    xc = array(ac[ac.columns[10:13]]).astype(float)
    x_calcium_mc = xc - (sum(xc)/len(xc))
    yc = array(bc[bc.columns[10:13]]).astype(float)
    y_calcium_mc = yc - (sum(yc)/len(yc))

    R, _ = align.rotation_matrix(x_mean_centered, y_mean_centered)

    # Rotate calciums based on the Rotational matrix from the previous alignment
    y_calcium_aligned = dot(y_calcium_mc,R)

    ##############Calciums##############
    # sort the coordinates to a point for ordering pairs
    x_sorted = array([list(i) for i in sorted(x_calcium_mc,key=lambda f: dist(f,x_mean_centered[0]))])
    y_sorted = array([list(i) for i in sorted(y_calcium_aligned,key=lambda f: dist(f,x_mean_centered[0]))])

    return round(rmsd(x_sorted,y_sorted),4)

import plotly
import plotly.graph_objs as go

def plot3d_ab(a,b):
    #colorlist = ['rgba(0,0,255,.1)' if i<start or i >stop else 'rgba(0,0,255,1)' for i in range(0,541*len(atomlist))]
    #colorlist2 = ['rgba(255,0,0,.1)' if i<start or i >stop else 'rgba(255,0,0,1)' for i in range(0,541*len(atomlist))]
    if len(a.columns)==3 and len(b.columns)==3:
        a = array(a[a.columns[10:13]]).astype(float)
        b = array(b[b.columns[10:13]]).astype(float)


    x,y,z = a[:,0],a[:,1],a[:,2]
    trace = go.Scatter3d(x=x,y=y,z=z,mode='lines+markers',
                         marker={'size': 5,"color":"blue"},line={'width':3},name="ptm structure")
    x,y,z = b[:,0],b[:,1],b[:,2]
    trace2 = go.Scatter3d(x=x,y=y,z=z,mode='lines+markers',
                          marker={'size': 5,"color":"red"},line={'width':3},name="noptm structure")
    
    # Configure the layout.
    layout = go.Layout(margin={'l': 0, 'r': 0, 'b': 0, 't': 0})
    data = [trace,trace2]
    plot_figure = go.Figure(data=data, layout=layout)
    # Render the plot.
    
    plot_figure.update_layout(autosize=True) # remove height=800
    #plot_figure.show(renderer="browser")  # remove display(fig)
    #plot_figure.write_html("plots/plot.html")
    # Adding vertical lines at z = 0
    
    plotly.offline.iplot(plot_figure)

def plot3d_a(a):
    #colorlist = ['rgba(0,0,255,.1)' if i<start or i >stop else 'rgba(0,0,255,1)' for i in range(0,541*len(atomlist))]
    #colorlist2 = ['rgba(255,0,0,.1)' if i<start or i >stop else 'rgba(255,0,0,1)' for i in range(0,541*len(atomlist))]
    
    a = array(a[a.columns[10:13]]).astype(float)

    x,y,z = a[:,0],a[:,1],a[:,2]
    trace = go.Scatter3d(x=x,y=y,z=z,mode='lines+markers',
                         marker={'size': 5,"color":"blue"},line={'width':3},name="ptm structure")
    
    
    # Configure the layout.
    layout = go.Layout(margin={'l': 0, 'r': 0, 'b': 0, 't': 0})
    data = [trace]
    plot_figure = go.Figure(data=data, layout=layout)
    # Render the plot.
    
    plot_figure.update_layout(autosize=True) # remove height=800
    #plot_figure.show(renderer="browser")  # remove display(fig)
    #plot_figure.write_html("plots/plot.html")
    # Adding vertical lines at z = 0
    
    plotly.offline.iplot(plot_figure)