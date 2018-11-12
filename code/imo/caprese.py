from rpy2.robjects.packages import importr
import pandas as pd
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
import numpy as np

T = importr('TRONCO')
class CAPRESE(object):
    def __init__(self):
        self.data = None
        self.driver_gens = None
        self.data_claen = None
        self.model = None


    def fit(self, X, driver_gens = None):
        r_dataframe = pandas2ri.py2ri(X)
        self.data = T.import_MAF(r_dataframe)
        if driver_gens is not None:
            self.driver_gens = r('c(' + str(driver_gens)[1:-1] + ')' )
            self.data_clean = T.events_selection(self.data, filter_in_names=self.driver_gens)
        else:
            self.data_claen = self.data

        self.model = T.tronco_caprese(self.data_clean)

    def draw(self, file_name=r('NA')):

        T.tronco_plot(self.model,
                      fontsize=12,
                      scale_nodes=0.6,
                      legend = False,
                      file=file_name)



