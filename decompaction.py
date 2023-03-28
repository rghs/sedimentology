# -*- coding: utf-8 -*-
"""
Created on Sun Nov 27 15:07:29 2022

@author: rghs

Functions to perform backstripping analysis. If you're looking at this, these
functions are probably still broken.
"""

import scipy.optimize as so
import numpy as np

def findY2(y2, ys, phic, c):
    y2 = ys + phic/c * (1 - np.e**(-c*y2))
    return y2