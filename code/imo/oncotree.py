import numpy as np
import pandas as pd
import networkx as nx


class OncoTree(object):
    """
    The tree model for oncogenesis by Desper, Richard, et al.

    This implementation is based on:

        Desper, R., Jiang, F., Kallioniemi, O. P., Moch, H., Papadimitriou, C. H., & Schäffer, A. A. (1999).
        Inferring tree models for oncogenesis from comparative genome hybridization data.
        Journal of computational biology, 6(1), 37-51. URL: https://www.ncbi.nlm.nih.gov/pubmed/10223663

    Parameters
    ----------
    min_mutations : integer or None, optional (default=None)
        Minimum number of mutations to consider. The genes which has mutated below this number
        will be dropped. If ``None``, no mutation is skipped.
    """

    def __init__(self, min_mutations=None):
        self.G = None
        self.edge_list = None
        self.optimum_branching = None
        self.min_mutations = min_mutations

    def create_init_graph(self, df):
        # for each gene, find the set of patients whom this gene has mutated in them and its size of this set
        cases = pd.DataFrame(df.groupby('Hugo_Symbol')['case_id'].apply(set))
        cases['count'] = cases['case_id'].apply(len)

        # drop low mutated genes if desired
        if self.min_mutations is not None:
            cases = cases[cases['count'] >= self.min_mutations]

        # for each pair of genes, find the size of the set of patients whom these two genes has both mutated in them
        cases['key'] = 0
        cases['symbol'] = cases.index
        self.edge_list = cases[['case_id', 'symbol', 'key']].merge(cases[['case_id', 'symbol', 'key']],
                                                                   how='outer', on='key')
        self.edge_list['count'] = self.edge_list.apply(lambda row: len(row['case_id_x'] & row['case_id_y'])
                                                       , axis=1)
        self.edge_list.drop(['case_id_x', 'case_id_y', 'key'], axis=1, inplace=True)

        n_cases = len(df['case_id'].unique())

        self.edge_list['weight'] = self.edge_list.apply(
            lambda row: np.log(row['count'] / n_cases) -
                        np.log(cases.loc[row['symbol_x']]['count'] / n_cases +
                               cases.loc[row['symbol_y']]['count'] / n_cases) -
                        np.log(cases.loc[row['symbol_y']]['count'] / n_cases)
            , axis=1)

        self.G = nx.from_pandas_edgelist(self.edge_list, source='symbol_x', target='symbol_y', edge_attr='weight',
                                         create_using=nx.DiGraph())
        return self.G

    def find_branching(self):
        self.optimum_branching = nx.algorithms.tree.branchings.maximum_spanning_arborescence(self.G)
        return self.optimum_branching

    def draw(self):
        pos = nx.drawing.nx_pydot.pydot_layout(self.optimum_branching, prog='dot')
        nx.draw(self.optimum_branching, pos, with_labels=True)
