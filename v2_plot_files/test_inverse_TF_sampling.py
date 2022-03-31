import numpy as np
import math
import random
import matplotlib.pyplot as plt

PI = math.pi

def sine_quantile_fun(p):
    return 1/PI*np.arccos(1-2*p)


samples = np.random.rand(int(1e6))
samples = sine_quantile_fun(samples)

plt.hist( samples, bins=100 )
plt.savefig( "figures/inverse_samplig.png" )
plt.close()
