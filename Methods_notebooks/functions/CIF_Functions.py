# Functions used in panda notebook and abstracted for use

# import modules
import pandas as pd
from scipy.stats import variation


def cif_line_to_list(line:str) -> list:
    """
    Function to seperate values by undisclosed amount of spaces
    Args:
        line (str): CIF style line containing multiple spaces

    Returns:
        list of strings of per column values
    """
    # split line by spaces, and remove (filter) all empty listitems
    # Also [:-1] to remove the \n at the end of every line 
    ls = list(filter(None,line.split(sep = " ")))[:-1]
    return ls


def cif_line_to_list_old(line:str) -> list:
    """
    Function to seperate values by undisclosed amount of spaces
    Args:
        line (str): CIF style line containing multiple spaces

    Returns:
        list of strings of per column values
    """
    
    line_list = []
    value = ""

    for char in line:
        if char != ' ':
            value+=char
        else:
            if len(value)== 0 or value == " ":
                pass
            else:
                line_list.append(value)
                value = ""

    return line_list


def make_df(filepath: str) -> tuple[pd.DataFrame,list[str]]:
    """Takes in a .cif file and returns dataframe object and residue sequence

    Args:
        filepath (str): 

    Returns:
        (pd.DataFrame, List[str] of residues)as tuple
    """
    # Empty lists to append to
    atom_list   = []
    columns     = []
    sequence    = []

    # Open the file
    with open(filepath, encoding="UTF-8") as file:
        for ln in file:
            # Find sequence lines
            if ln.startswith("1 n"):
                sequence.append(cif_line_to_list(ln)[2])
            # Find atom lines
            elif ln.startswith("ATOM"):
                atom_list.append(cif_line_to_list(ln))
            # Get column names
            elif ln.startswith("_atom_site."):
                columns.append(ln[11:-1])

    # Make the pandas Dataframe
    df = pd.DataFrame(atom_list)
    df.columns= columns

    return df,sequence


def make_df_list(filepath: str) -> tuple[list[pd.DataFrame],list[list[str]]]:
    """Takes in a .cif file and returns list of dataframe objects 
    and residue sequences per found protein
    
    Args:
        filepath (str): 

    Returns:
        list of (pd.DataFrame, List[str] of residues) as tuple
    """
    # Getting the last line to check N_proteins
    with open(filepath, encoding="UTF-8") as file:
        atomline = ""
        for line in file:
            if line.startswith("ATOM"):
                atomline = line

        last_line = atomline 
    
    N_proteins = ord(cif_line_to_list(last_line)[-2])-64

    
    # Empty lists to append to
    atom_list   = []
    columns     = []
    sequences    = [[] for _ in range(N_proteins)]
    model_list = [f"{i} n" for i in range(1,N_proteins+1)]

    # Open the file
    with open(filepath, encoding="UTF-8") as file:
        for ln in file:
            # Find sequence lines
            if ln[0:3] in model_list:
                model = int(ln[0])
                sequences[model-1].append(cif_line_to_list(ln)[2])
            # Find atom lines
            elif ln.startswith("ATOM"):
                atom_list.append(cif_line_to_list(ln))
            # Get column names
            elif ln.startswith("_atom_site."):
                columns.append(ln[11:-1])

    # Make the pandas Dataframe
    df = pd.DataFrame(atom_list)
    df.columns= columns

    df_list = []
    # The rows are still all models in the ATOM list so we seperate those by Model
    for i in range(1,N_proteins+1):
        df_list.append(df[(df["label_asym_id"] == chr(i+64))])

    return df_list,sequences


def group_mean_variation(tpl:tuple)-> pd.DataFrame:
    '''Takes a dataframe of ATOM scores and groups them by residue
    Args:
        tuple of a dataframe and its sequence

    Returns:
        list of strings of per column values
    
    '''
    #unpack tuple
    df, sequence = tpl

    # Group residue numbers (atoms belonging to the same residue)
    id_score_df = df[["label_seq_id","B_iso_or_equiv"]].astype(float)
    grouped = id_score_df.groupby("label_seq_id")

    # Show the mean of those groups
    sp_residue = grouped.mean()

    # Show various statistics on these groups
    sp_residue["Residue_name"]  = sequence
    sp_residue["Variation"]     = grouped.agg(lambda x: variation(x))
    sp_residue["n_Atoms"]       = grouped.agg(lambda x: len(x))
    #sp_residue["variance"]     = grouped.var()
    #sp_residue["stdev"]        = grouped.std()
    return sp_residue # type: ignore


def make_df_xr(filepath: str) -> tuple[list[pd.DataFrame],list[list[str]]]:
    """Takes in a .cif file and returns dataframe object and residue sequence

    Args:
        filepath (str): 

    Returns:
        (pd.DataFrame, List[str] of residues)as tuple
    """
    # Getting the last line to check Number of proteins (chains)
    with open(filepath, encoding="UTF-8") as file:
        atomline = ""
        for line in file:
            if line.startswith("ATOM"):
                atomline = line
            pass

        last_line = atomline
    N_proteins = ord(cif_line_to_list(last_line)[6])-64

    # Empty lists to append to
    atom_list   = []
    columns     = []
    sequences    = [[] for _ in range(N_proteins)]
    model_list = [f"{i} n" for i in range(1,N_proteins+1)]

    # Open the file
    with open(filepath, encoding="UTF-8") as file:
        for ln in file:
            # Find sequence lines
            if ln[0:3] in model_list:
                model = int(ln[0])
                sequences[model-1].append(cif_line_to_list(ln)[2])
            # Find atom lines
            elif ln.startswith("ATOM"):
                atom_list.append(cif_line_to_list(ln))
            # Get column names
            elif ln.startswith("_atom_site."):
                columns.append(ln[11:-1])

    # Make the pandas Dataframe
    df = pd.DataFrame(atom_list)
    df.columns= columns
    #print(columns)
    df_list = []

    # The rows are still all models in the ATOM list so we seperate those by Model
    for i in range(1,N_proteins+1):
        df_list.append(df[(df["label_asym_id "] == chr(i+64))])

    return df_list,sequences


def make_df_new(filepath: str) -> tuple[list[pd.DataFrame],list[list[str]]]:
    """Takes in a .cif file and returns dataframe object and residue sequence per chain found (ABCD)

    Args:
        filepath (str): 

    Returns:
        (pd.DataFrame, List[str] of residues)as tuple
    """
    # Getting the last line to check Number of proteins (chains)
    with open(filepath, encoding="UTF-8") as file:
        atomline = ""
        for line in file:
            if line.startswith("ATOM"):
                atomline = line

        last_line = atomline
    
    # Setting number to ordinal of last line (A=1,B=2,etc)
    N_proteins = ord(cif_line_to_list(last_line)[6])-64

    # Empty lists to append to
    atom_list   = []
    columns     = []
    
    # Open the file
    with open(filepath, encoding="UTF-8") as file:
        for ln in file:
            # Find atom lines
            if ln.startswith("ATOM"):
                atom_list.append(cif_line_to_list(ln))
            # Get column names
            elif ln.startswith("_atom_site."):
                columns.append(ln[11:-1])

    # Make the pandas Dataframe
    df = pd.DataFrame(atom_list)

    # Crystals downloaded somehow have a whitespace after the column name so il clip that off
    if str(columns[0][-1]) == ' ':
        columns = [column[:-1] for column in columns]
    df.columns= columns

    #empty list to separate different models in
    df_list = []
    sequences = []

    # The rows are still all models in the ATOM list so we seperate those by Model
    for i in range(1,N_proteins+1):
        # Dataframes seperated by model A,B,C etc
        df_to_add = df[(df[columns[6]] == chr(i+64))]
        # Take the residue column (3) from each C-Alpha to fetch the protein sequence
        sequences.append(list(df_to_add[df_to_add[columns[3]]=="CA"][columns[5]]))

        df_list.append(df_to_add)

    return df_list,sequences

def hetatm_df(filepath: str) -> pd.DataFrame:
    """Takes in a .cif file and returns dataframe object and residue sequence

    Args:
        filepath (str): 

    Returns:
        (pd.DataFrame, List[str] of residues)as tuple
    """
    # Empty lists to append to
    atom_list   = []
    columns     = []
   
    # Open the file
    with open(filepath, encoding="UTF-8") as file:
        for ln in file:
            # Find sequence lines
            if ln.startswith("HETATM"):
                atom_list.append(cif_line_to_list(ln))
            # Get column names
            elif ln.startswith("_atom_site."):
                columns.append(ln[11:-1])

    # Make the pandas Dataframe
    df = pd.DataFrame(atom_list)
    
    # Crystals downloaded somehow have a whitespace after the column name so il clip that off
    if str(columns[0][-1]) == ' ':
        columns = [column[:-1] for column in columns]
    df.columns= columns

    return df

def make_df_atomlist(path,atom_list):
    raw_df = make_df_new(path)[0][0]
    im_df = raw_df[raw_df["label_atom_id"].isin(atom_list)]
    return im_df[(im_df['label_alt_id'] == ".")|(im_df['label_alt_id']=="A")]


def make_df_full(filepath: str) -> pd.DataFrame:
    """Takes in a .cif file and returns dataframe object and residue sequence per chain found (ABCD)

    Args:
        filepath (str): 

    Returns:
        (pd.DataFrame, List[str] of residues)as tuple
    """

    # Empty lists to append to
    atom_list   = []
    columns     = []
    
    # Open the file
    with open(filepath, encoding="UTF-8") as file:
        for ln in file:
            # Find atom lines
            if ln.startswith("ATOM") or ln.startswith("HETATM"):
                atom_list.append(cif_line_to_list(ln))
            # Get column names
            elif ln.startswith("_atom_site."):
                columns.append(ln[11:-1])

    # Make the pandas Dataframe
    df = pd.DataFrame(atom_list)

    # Crystals downloaded somehow have a whitespace after the column name so ill clip that off
    if str(columns[0][-1]) == ' ':
        columns = [column[:-1] for column in columns]
    df.columns= columns
    
    # We select only the A alternative id
    df = df[(df['label_alt_id'] == ".")|(df['label_alt_id']=="A")]

    # We select only the given chain "A"
    
    
    return df