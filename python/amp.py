#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from pylab import *


def prior(r, s, t):
    rho, m_pr, v_pr = t
    
    m_eff = (m_pr * s + r * v_pr) / (s + v_pr)
    v_eff = v_pr * s / (s + v_pr)
    z = ((1. - rho) / rho) * sqrt(v_pr / v_eff) * \
       exp(-.5 * ( r ** 2 / s - (m_pr - r) ** 2 / (s + v_pr) ))
    
    a = m_eff / (z + 1)
    c = z * a ** 2 + v_eff / (z + 1)
    
    return a, c

def run_amp(y, F, x, d, prior_prmts, t_max):
    m, n = F.shape
    sqrF = F * F
    
    a, c = zeros(n), ones(n)
    w, v = y, ones(m)
    
    mse = []
    for t in range(t_max):
        w = F.dot(a) - sqrF.dot(c) * (y - w) / (d + v)
        v = sqrF.dot(c)
        s = 1. / sqrF.T.dot( 1. / (d + v) )
        r = a + s * F.T.dot( (y - w) / (d + v) )
        
        a_old = a
        a, c = prior(r, s, prior_prmts)
        diff = mean(abs(a - a_old))

        # Uncomment this line to learn \Delta iteratively.
        #d = d * sum(((y - w) / (d + v)) ** 2) / sum(1. / (d + v))

        mse += [mean((a - x) ** 2)]
    mse = array(mse)
    
    return a, c, mse, diff

def run_swamp(y, F, x, d, prior_prmts, t_max):
    m, n = F.shape
    sqrF = F * F
    
    a, c = zeros(n), ones(n)
    r, s = zeros(n), zeros(n)
    w, v = y, ones(m)
    
    mse = []
    for t in range(t_max):
        diff = 0

        g = (y - w) / (d + v)
        v = sqrF.dot(c)
        w = F.dot(a) - v * g
        for i in permutation(n):
            s[i] = 1. / sqrF[:, i].dot( 1. / (d + v) )
            r[i] = a[i] + s[i] * F[:, i].dot( (y - w) / (d + v) )
        
            a_old, c_old = a[i], c[i]
            a[i], c[i] = prior(r[i], s[i], prior_prmts)
            diff += (a[i] - a_old)

            v_old = v;
            v += sqrF[:, i].dot(c[i] - c_old)
            w += F[:, i].dot(a[i] - a_old) - (v - v_old) * g
        # We are learning \Delta iteratively here; we don't have to, though!
        d = d * sum(((y - w) / (d + v)) ** 2) / sum(1. / (d + v))

        mse += [mean((a - x) ** 2)]
    mse = array(mse)
    
    return a, c, mse, diff
