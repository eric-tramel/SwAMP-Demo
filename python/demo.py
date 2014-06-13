#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from pylab import *
import time

from amp import run_amp, run_swamp

seterr(all = 'ignore')


def sample_instance(n, m, k, delta, gamma):
    x = zeros(n)
    supp = np.random.choice(n, k, replace = False)
    x[supp] = normal(size = k)
    
    F = float(gamma) / n + normal(size = (m, n)) / sqrt(n);
    y = normal(F.dot(x), sqrt(delta))
    
    return x, F, y

def nzm(gamma):
    # Instance's parameters
    n = 1024
    alpha = 0.72
    rho, m_pr, v_pr = 0.44, 0.0, 1.0
    delta = 1e-8

    print r' * Parameters: N = %d, \rho = %.2f, \alpha = %.2f, \Delta = %.2g, \gamma = %d' % \
            (n, rho, alpha, delta, gamma)

    prior_prmts = array([rho, m_pr, v_pr])

    # Pre-processing
    m = ceil(alpha * n)
    k = ceil(prior_prmts[0] * n)

    # Gen. instance
    x, F, y = sample_instance(n, m, k, delta, gamma)

    # Algorithm's parameters
    t_max = 200

    # Run algorithm
    print ' - Running SwAMP...'
    t = time.time()
    a_sw, c_sw, mse_sw, diff_sw = run_swamp(y, F, x, 1.0, prior_prmts, t_max)
    elapsed = time.time() - t
    print r'  Elapsed time: %.2fs, final MSE: %.4g' % (elapsed, mse_sw[-1])

    print ' - Running AMP...'
    t = time.time()
    a_amp, c_amp, mse_amp, diff_amp = run_amp(y, F, x, delta, prior_prmts, t_max)
    elapsed = time.time() - t
    print r'  Elapsed time: %.2fs, final MSE: %.4g' % (elapsed, mse_amp[-1])

    # Plot results
    subplot(2, 1, 1)
    plot(x, 'ko', ms = 7, mfc = 'none'); plot(a_sw, 'rx', ms = 5); plot(a_amp, 'b+', ms = 5)
    xlim([0, len(x)]); xlabel('i'); ylabel('x')

    subplot(2, 2, 3)
    semilogy(mse_sw, 'r')
    xlim([0, len(mse_sw)]); xlabel('SwAMP iter.'); ylabel('MSE')

    subplot(2, 2, 4)
    semilogy(mse_amp, 'r')
    xlim([0, len(mse_amp)]); xlabel('AMP iter.')

    gcf().set_tight_layout(True)
    show()

print ' Welcome to the SwAMP demo!'
print " WARNING: This version of the demo is 'conceptual': it's not at all"
print "   optimized, and presents only a single example. If you have MATLAB"
print "   installed on your computer, please use the other version."
raw_input(' (Press any key.)')
print
print " Let's run both AMP and SwAMP for i.i.d. Gaussian matrices with zero"
print "   mean and 1 / N variance."
print
nzm(0)
print
print " As expected, both work and reach the same MSE. We'll now try setting"
print "   the projector mean to 20 / N."
raw_input(' (Press any key.)')
print
nzm(20)
print
print " In this case, AMP diverges, while SwAMP remains converging to a good MSE."
