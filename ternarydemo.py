# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 12:22:56 2021

@author: RGHS
"""

import sys
sys.path.append('.')

import numpy as np
import graphs as gr
import pandas as pd

import plotly.io as pio

pio.renderers.default = 'svg'

data = pd.read_csv('D:/code/input/ternplotsamples.csv')

qtfl_qt = np.array((data.Qp+data.Qm)/(data.Qp+data.Qm+data.P+data.K+data.L))
qtfl_f = np.array((data.P+data.K)/(data.Qp+data.Qm+data.P+data.K+data.L))
qtfl_l = np.array((data.L)/(data.Qp+data.Qm+data.P+data.K+data.L))

qmflt_qm = np.array((data.Qm)/(data.Qp+data.Qm+data.P+data.K+data.L))
qmflt_f = np.array((data.P+data.K)/(data.Qp+data.Qm+data.P+data.K+data.L))
qmflt_lt = np.array((data.L+data.Qp)/(data.Qp+data.Qm+data.P+data.K+data.L))

data = data.join(pd.DataFrame({'qtfl_qt'    :qtfl_qt*100,
                               'qtfl_f'     :qtfl_f*100,
                               'qtfl_l'     :qtfl_l*100,
                               'qmflt_qm'   :qmflt_qm*100,
                               'qmflt_f'    :qmflt_f*100,
                               'qmflt_lt'   :qmflt_lt*100}))

ternary = gr.qtfl(data, 'qtfl_qt','qtfl_f','qtfl_l', label = 'sample')
ternary.show()