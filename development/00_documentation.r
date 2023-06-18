################################################################################
# cyclar
# by Peter Hellstr?m, Department of Zoology, Stockholm University
# Last major update: October 17th 2011
# Minor update: September 17th 2020
################################################################################
# The cyclar "package" contains:
# Functions for simulation and estimation of time series
# Mainly based on the AR(2)-process

# Different sets of functions
# 1) Investigate the dynamics of an AR(2)-process
	# a) Plot dynamics in (1 + a1) / a2 space
	# b) Functions for calculating (1 + a1) and a2, given k (period) and ipv (intrinsic process variance) - or vice versa.
# 2) Simulate a stationary AR(2)-process (uncorrelated and correlated)
# 3) Estimate parameters

# Various time series functions (summary statistics & detrending)
# Autocorrelation functions
# Markov chains (simulate and phase/state transitions)
# Various functions
# Time-varying parameters

################################################################################
# KNOWN ISSUES
# a1 in this code (mostly) refers to (1 + a1), written here as a1 for simplicity!
# But beware of the difference between a1 and (1 + a1) when dealing with Rt and xt!

# Functions with "phi"-arguments
# ar2.ipv2
# ar2.acf2
# ar2.period2
# ar2.plot.simple2
# ar2.arrows2
# ar2.intercept2
# ar2.k2
# ar2.parabola2
# ts.spec.ar22
