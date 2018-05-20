import pandas as pd


def read_maf(filename, skip_first_lines=True):
    """
    Mutation Annotation Format (MAF) parser

    Parameters
    ----------
    filename : string
        The name of the tab-delimited text (.tsv) file containing the mutation information
    skip_first_lines : boolean (default=True)
        Whether to skip the first 5 lines which are the file information

    Returns
    -------
    pandas.DataFrame
        The data frame containing the result of parsing
    """
    df = pd.read_table(filename, skiprows=5 if skip_first_lines else None)
    return df
