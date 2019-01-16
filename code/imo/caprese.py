from rpy2.robjects.packages import importr
import pandas as pd
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
import numpy as np
from rpy2.interactive import process_revents
from rpy2.robjects.packages import importr
process_revents.start()
import rpy2
import IPython
from rpy2.robjects.lib import grdevices
import IPython

graphics = importr("graphics")

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

        with rpy2.robjects.lib.grdevices.render_to_bytesio(grdevices.png, width=1024, height=896, res=150) as img:
            T.tronco_plot(self.model,
                          fontsize=12,
                          scale_nodes=0.6,
                          legend = False,
                          file=file_name)
        IPython.display.display(IPython.display.Image(data=img.getvalue(), format='png', embed=True))



