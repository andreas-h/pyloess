#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pickle

import numpy
from scipy.stats.mstats import kendalltau_seasonal

# set up input parameters

filename_ts = "/home2/hilboll/data/pickle/20100920_trends_levelshift_cityzen/Pearl River Delta_seas-const.pkl"
with open(filename_ts, 'rb') as _fd:
    ts = pickle.load(_fd)

prd = ts[2][0].fill_missing_dates()

y = numpy.ma.masked_all(prd.shape[0] + 3)
y[3:] = prd
#
np = 12         # period of seasonal component
ns = 7          # length of seasonal smoother
nt = None       # length of trend smoother
nl = None       # length of low-pass filter
isdeg = 1       # Degree of locally-fitted polynomial in seasonal smoothing.
itdeg = 1       # Degree of locally-fitted polynomial in trend smoothing.
ildeg = 1       # Degree of locally-fitted polynomial in low-pass smoothing.
nsjump = None   # Skipping value for seasonal smoothing.
ntjump = 1      # Skipping value for trend smoothing. If None, ntjump= 0.1*nt
nljump = 1      # Skipping value for low-pass smoothing. If None, nljump= 0.1*nl
robust = True   # Flag indicating whether robust fitting should be performed.
ni = None       # Number of loops for updating the seasonal and trend  components.
no = 0          # Number of iterations of robust fitting. The value of no should

def run_ndarray():
    from pyloess import stl
    y = prd.compressed()
    res = stl(y, np, ns, nt, nl, isdeg, itdeg, ildeg, nsjump, ntjump, nljump, robust, ni, no)

    figure()
    plot(y, "k")
    plot(res.trend, "b")
    plot(res.seasonal, "g")
    plot(res.residuals, "r")
    figure()
    plot(res.weights)


def run_maskedarray():
    from pyloess.mpyloess import stl
    res = stl(y)#, np=np, ns=ns, nt=nt, nl=nl, isdeg=isdeg, itdeg=itdeg, ildeg=ildeg, nsjump=nsjump, ntjump=ntjump, nljump=nljump, robust=robust, ni=ni, no=no)

    figure()
    plot(y, "k")
    plot(res.trend, "b")
    plot(res.seasonal, "g")
    plot(res.residuals, "r")
    figure()
    plot(res.weights)

