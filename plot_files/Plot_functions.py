import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import csv
import matplotlib as mpl


class animate_density:
    def __init__(self, ax, positions ):
        self.ax = ax
        self.positions = positions
        
