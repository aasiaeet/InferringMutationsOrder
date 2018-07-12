import networkx as nx
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

    @classmethod
    def __get_all_genes(cls, u):
        return u.index

    def _reduce_events(self, u):
        pass

    def _find_structure(self, u):
        self.G = nx.DiGraph()
        for f in self.events:
            for e in self.events:
                if u.apply(lambda g: g & {e, f} == {f}).sum() <= self.epsilon:
                    self.G.add_edge(e, f)
        return self.G

    def _find_probabilities(self, u):
        theta = {}
        for e in self.G.nodes:
            below = u.apply(lambda g: set(self.G.predecessors(e)).issubset(g)).sum()
            theta[e] = round(u.apply(lambda g: e in g).sum() / below, 2)
        nx.set_node_attributes(self.G, theta, 'theta')
        return theta

    def fit(self, X):
        u = CBN.__df_to_prob_dist(X)

        if self.events is None:
            self.events = CBN.__get_all_genes(u)
        if 0 < self.epsilon < 1:
            self.epsilon = int(u.shape[0] * self.epsilon)

        self._reduce_events(u)
        self._find_structure(u)
        self._find_probabilities(u)

    def draw(self, figsize=(20, 15), with_node_probs=True):
        plt.figure(figsize=figsize)
        pos = nx.drawing.nx_pydot.pydot_layout(self.G, prog='dot')
        labels = {}
        if with_node_probs:
            for e in self.G.nodes:
                labels[e] = e + ':' + str(nx.get_node_attributes(self.G, 'theta')[e])
        nx.draw(self.G, pos, labels=labels)

