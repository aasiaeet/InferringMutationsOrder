from rpy2.robjects.packages import importr
import pandas as pd
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
import numpy as np

T = importr('TRONCO')
class CAPRI(object):
    def __init__(self, nboot=5, boot_seed=12345):
        self.data = None
        self.driver_gens = None
        self.data_claen = None
        self.model = None
        self.boot_seed = boot_seed
        self.nboot = nboot



    def fit(self, X, driver_gens = None):
        r_dataframe = pandas2ri.py2ri(X)
        self.data = T.import_MAF(r_dataframe)
        if driver_gens is not None:
            self.driver_gens = r('c(' + str(driver_gens)[1:-1] + ')' )
            self.data_clean = T.events_selection(self.data, filter_in_names=self.driver_gens)
        else:
            self.data_claen = self.data

        self.model = T.tronco_capri(self.data_clean, boot_seed=self.boot_seed, nboot=self.nboot)


    def draw(self,file_name = r('NA')):

        T.tronco_plot(self.model,
                      fontsize=12,
                      scale_nodes=0.6,
                      confidence=r("c('tp', 'pr', 'hg')"),
                      height_logic=0.25,
                      legend_cex=0.35,
                      pathways=list(self.driver_gens),
                      label_edge_size=10,
                      file = file_name)


