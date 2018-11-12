import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


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

        # compute the number of cases
        n_cases = len(df['case_id'].unique())

        # compute the edge weights of the graph using the formula:
        # w_ij = log(p_ij) - log(p_i + p_j) - log(p_j)
        self.edge_list['weight'] = self.edge_list.apply(
            lambda row: np.log(row['count'] / n_cases) -
                        np.log(cases.loc[row['symbol_x']]['count'] / n_cases +
                               cases.loc[row['symbol_y']]['count'] / n_cases) -
                        np.log(cases.loc[row['symbol_y']]['count'] / n_cases)
            , axis=1)

        # create a networkx graph from the edgelist we just created
        self.G = nx.from_pandas_edgelist(self.edge_list, source='symbol_x', target='symbol_y', edge_attr='weight',
                                         create_using=nx.DiGraph())
        return self.G

    def find_branching(self):
        """
        Find the oncogenetic graph using Edmond's branching algorithm. Needs ``create_init_graph`` function
        to be called before using this function.

        Returns
        -------
        branching: networkx.DiGraph
            the oncogenetic tree

        """
        self.optimum_branching = nx.algorithms.tree.branchings.maximum_spanning_arborescence(self.G)
        return self.optimum_branching

    def fit(self, df):
        """
        Runs the Desper's algorithm on the provided data frame.

        Parameters
        ----------
        df: pandas.DataFrame
            The data matrix. Should have ``Hugo_Symbol`` and ``case_id`` columns.
        """
        self.create_init_graph(df)
        self.find_branching()

    def draw(self,
             figsize=(20, 15),
             with_edges=True,
             node_color= 'white',
             node_size= 10000,
             font_color='black',
             font_size = 20,
             width_ = 3,
             arrowsize_ = 30,
             edge_label_size = 20,
             edge_node_color = 'black'
             ):  # TODO move this function to utils
        """
        Draws the oncogenetic tree. Needs ``find_branching`` to be called before drawing.
        """
        # set the figure size
        plt.figure(figsize=figsize)

        # save the edge weights for later use
        edge_weights = {(u, v): round(self.optimum_branching[u][v]['weight'], 2) for u, v in
                        self.optimum_branching.edges}
        # change all edge weights of the branching to 1 to get a nice hierarchical tree in drawing
        for u, v in self.optimum_branching.edges:
            self.optimum_branching[u][v]['weight'] = 1
        pos = nx.drawing.nx_pydot.pydot_layout(self.optimum_branching, prog='dot')
        nx.draw(self.optimum_branching,
                pos,
                with_labels = True,
                node_color = node_color,
                node_size = node_size,
                font_color = font_color,
                font_size = font_size,
                width = width_,
                arrowsize = arrowsize_,
                edgecolors = edge_node_color)

        if with_edges:
            nx.draw_networkx_edge_labels(self.optimum_branching,
                                         pos,
                                         edge_labels=edge_weights,
                                         font_size = edge_label_size)

        # reset the edge weights to their original weights
        for u, v in edge_weights:
            self.optimum_branching[u][v]['weight'] = edge_weights[(u, v)]
