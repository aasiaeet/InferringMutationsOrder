from rpy2.robjects.packages import importr
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
from rpy2.robjects.lib import grdevices
import rpy2
import pandas as pd
Rplot = rpy2.robjects.r('plot')
import IPython
T = importr('Rtreemix')
importr('graph')
importr('Rgraphviz')
importr('grid')

Rplot = rpy2.robjects.r('plot')
class Treemix(object):
    def __init__(self, nboot=5, boot_seed=12345):
        self.data = None
        self.driver_gens = None
        self.data_claen = None
        self.model = None


    def _read_data(self, X, drivers):
        X = X[X.Hugo_Symbol.isin(drivers)]
        gr = pd.DataFrame(X.groupby('case_id')['Hugo_Symbol'].apply(list))
        data = dict()
        for i in drivers:
            person = dict()
            for j in gr.index:
                person[j] = i in list(gr[gr.index == j]['Hugo_Symbol'])[0]
            data[i] = person

        data = pd.DataFrame(data)
        data2 = data.apply(lambda x: x.astype('int'))
        data2.index = range(len(data2.index))
        cbind = rpy2.robjects.r('cbind')
        CR = rpy2.robjects.r('c')
        all_pataint = cbind()
        for i in data2.values.T.tolist():
            j = list(i)
            all_pataint = cbind(all_pataint, CR(*j))

        clean_data = rpy2.robjects.r('new')("RtreemixData", Sample=all_pataint)

        for i in range(len(drivers) + 1):
            T.Events(clean_data)[i] = CR('zero', CR(list(data2)))[i]
        return clean_data
    def fit(self, X, driver_gens = None):
        self.data_claen = self._read_data(X=X, drivers=driver_gens)
        self.model = model = T.fit(data = self.data_claen, K = 2, noise = True)




    def draw(self):
        with rpy2.robjects.lib.grdevices.render_to_bytesio(grdevices.png, width=1024, height=896, res=150) as img:
            Rplot(self.model, fontSize=14)

        IPython.display.display(IPython.display.Image(data=img.getvalue(), format='png', embed=True))


