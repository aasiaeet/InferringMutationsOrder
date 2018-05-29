import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
class OncoTreeWithFpFn(object):
    """
    The tree model for oncogenesis by Szabo, Aniko, et al.

    This implementation is based on:

        Beerenwinkel, N., Schwarz, R. F., Gerstung, M., & Markowetz, F. (2015).
        Cancer evolution: Mathematical models and computational inference.
        Systematic Biology, 64(1), e1–e25. https://doi.org/10.1093/sysbio/syu081

    Parameters
    ----------
    min_mutations : integer or None, optional (default=None)
        Minimum number of mutations to consider. The genes which has mutated below this number
        will be dropped. If ``None``, no mutation is skipped.
    """

    def __init__(self, min_mutations=None):
        self.edge_list = None
        self.opontimum_branching = None
        self.min_mutations = min_mutations
        self.mutation_freq = None
        self.branching_edge = None

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

        # compute the number of cases
        n_cases = len(df['case_id'].unique())
        self.edge_list = self.edge_list[self.edge_list['symbol_x'] != self.edge_list['symbol_y']]

        # compute the edge weights of the graph using the formula:
        # w_ij = log(p_ij) - log(p_i + p_j) - log(p_j)
        self.edge_list['weight'] = self.edge_list.apply(
            lambda row: np.log(row['count'] / n_cases) -
                        np.log(cases.loc[row['symbol_x']]['count'] / n_cases +
                               cases.loc[row['symbol_y']]['count'] / n_cases) -
                        np.log(cases.loc[row['symbol_y']]['count'] / n_cases)
            , axis=1)
        # a list of mutation sorted by frequency
        self.mutation_freq = list(cases.sort_values('count')['symbol'])

    def find_branching(self):

        """

        find the minimum branching tree by finding the parent of every nodes


        :return: optimum branching
        """
        tree = list()
        all_edge = self.edge_list.copy()
        all_edge.sort_values('weight', ascending=False, inplace=True)
        for i in self.mutation_freq:
            try:
                tree.append(all_edge[all_edge.symbol_y == i].iloc[0])
                all_edge = all_edge[all_edge.symbol_x != i]
            except:
                print("Warning: there is a mutation that has no parent.")
        self.branching_edge = pd.DataFrame(tree)
        self.optimum_branching = nx.from_pandas_edgelist(self.branching_edge, source='symbol_x', target='symbol_y',
                            edge_attr='weight', create_using=nx.DiGraph())

        return self.optimum_branching


    def fit(self, df):
        """
        Runs algorithm on the provided data frame.


        :param df: pandas.DataFrame
            The data matrix. Should have ``Hugo_Symbol`` and ``case_id`` columns.

        """
        self.create_init_graph(df)
        self.find_branching()


    def draw(self, figsize= (30,20), weight = False):
        """
        Plot the graph

        :param figsize:
        :param weight: show the weight of edge
        """
        plt.figure(figsize = figsize)
        if weight == False:
            pos = nx.drawing.nx_pydot.pydot_layout(self.optimum_branching, prog='dot')
            nx.draw(self.optimum_branching, pos, with_labels=True)
        else:
            pos = nx.drawing.nx_pydot.pydot_layout(self.optimum_branching, prog='dot')
            nx.draw_networkx_edge_labels(self.optimum_branching, pos,
                                         edge_labels={(u, v): round(self.optimum_branching.get_edge_data(u, v)['weight'], 2) for u, v in
                                         self.optimum_branching.edges});

            nx.draw(self.optimum_branching, pos, with_labels=True)
