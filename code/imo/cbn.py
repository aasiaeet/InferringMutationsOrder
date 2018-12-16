import networkx as nx
import numpy as np
from matplotlib import pyplot as plt


class CBN(object):

    def __init__(self, epsilon=None, genes=None):
        self.epsilon = epsilon
        if epsilon is None:
            self.epsilon = 0
        self.events = genes
        self.G = None

    @classmethod
    def __df_to_prob_dist(cls, df):
        return df.groupby('case_id')['Hugo_Symbol'].apply(set)

    def _reduce_events(self, u):
        pass

    def _find_structure(self, u):
        self.G = nx.DiGraph()
        for f in self.events:
            for e in self.events:
                n_errors = u.apply(lambda g: g & {e, f} == {f}).sum()
                if n_errors <= self.epsilon:
                    self.G.add_edge(e, f)
        return self.G

    def _find_probabilities(self, u):
        theta = {}
        nodes = np.copy(self.G.nodes)
        for e in nodes:
            below = u.apply(lambda g: set(self.G.predecessors(e)).issubset(g)).sum()
            if below == 0:
                self.G.remove_node(e)
            else:
                theta[e] = round(u.apply(lambda g: e in g).sum() / below, 3)
        nx.set_node_attributes(self.G, theta, 'theta')
        return theta

    def fit(self, X):
        all_genes = X['Hugo_Symbol'].unique()
        u = CBN.__df_to_prob_dist(X)

        if self.events is None:
            self.events = all_genes
        self.events = [e for e in self.events if e in all_genes]
        if 0 < self.epsilon < 1:
            self.epsilon = int(u.shape[0] * self.epsilon)

        self._reduce_events(u)
        self._find_structure(u)
        self._find_probabilities(u)

    def draw(self,
             figsize=(20, 15),
             with_node_probs=True,
             node_color='white',
             node_size=10000,
             font_color='black',
             font_size=20,
             width_=3,
             arrowsize_=30,
             edge_label_size=20,
             edge_node_color='black'):
        plt.figure(figsize=figsize)
        pos = nx.drawing.nx_pydot.pydot_layout(self.G, prog='dot')
        labels = {}
        if with_node_probs:
            for e in self.G.nodes:
                labels[e] = e + ':' + str(nx.get_node_attributes(self.G, 'theta')[e])
        nx.draw(self.G,
                pos,
                labels=labels,
                node_color=node_color,
                node_size=node_size,
                font_color=font_color,
                font_size=font_size,
                width=width_,
                arrowsize=arrowsize_,
                edgecolors=edge_node_color)
