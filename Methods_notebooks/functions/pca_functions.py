"PCA analysis condensed in a single function"

from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis import pca, align

import MDAnalysis as mda
import pandas as pd
import numpy as np


def pca_on_traj(list_of_lists,listnames,selec="index 0:1700"):
    """list of lists is list of pdb PATHS, first list is used for PCA 
    
    """
    universe_lists = [[mda.Universe(path).select_atoms(selec) for  path in file_list] for file_list in list_of_lists]
    ref_atoms = universe_lists[0][0]

    # Align all to the first frame and collect aligned coordinates
    coords_all =[]
    for universe_list in universe_lists:
        aligned_coords = []
        for mobile_atoms in universe_list:
            align.alignto(mobile_atoms, ref_atoms, select=selec)
            aligned_coords.append(mobile_atoms.positions.copy())
        coords_all.append(np.stack(aligned_coords, axis=0))

    # run pca on first list only
    u = mda.Universe(list_of_lists[0][0])
    u.trajectory = MemoryReader(coords_all[0])
    pc = pca.PCA(u, select=selec,align=True, mean=None,n_components=10).run()

    dfs = []
    for structures,coords,ls in zip(list_of_lists,coords_all,listnames):
        # project structure and its traj(coords) on pca
        u = mda.Universe(structures[0])
        u.trajectory = MemoryReader(coords)
        backbone = u.select_atoms(selec)
        transformed = pc.transform(backbone, n_components=3)

        df = pd.DataFrame(transformed,
                        columns=['PC{}'.format(i+1) for i in range(3)])
        df['Frame'] = df.index * u.trajectory.dt
        df["name"]= ['_'.join(path.split("fold_")[-1].split("_model")[0].split("_")[-2:]) for path in structures]
        df["list"]= ls
        dfs.append(df)

    return pd.concat(dfs)

