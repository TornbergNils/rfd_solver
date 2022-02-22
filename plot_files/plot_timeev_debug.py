import numpy as np
import matplotlib.pyplot as plt

time_ev = np.loadtxt( "./data/time_ev_at_pt.csv", delimiter=',' )

plt.plot(time_ev[:,1], time_ev[:,0] )
plt.plot(time_ev[:,1], np.sin(2.0* time_ev[:,1] + np.pi/1.9)*2.0, '--')
plt.show()
plt.close()
