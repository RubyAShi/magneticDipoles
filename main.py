from singleDipole import SquidLayout
from singleDipole import SingleDipole
from RandomDipoles import RandomDipoles
from singleDipole import twoAnglesConvolve
import numpy as np

# to generate 3um SQUID layout
#a = SquidLayout()
#b = a.IBM3um_mask()

# to scale a 3um SQUID layout according to aligning angle
#a = SquidLayout()
#b = a.IBM3um_mask()
#b = a.y_scales()

# to generate magnetic field of a number of randomly aligned dipoles
#c = RandomDipoles(Rx=40, Ry=40, num=100)
#c.generate_ran()
#d = c.B_z_ran()

# to convolve magnetic field from above with SQUID pickup loop
#a = SquidLayout()
#a.IBM3um_mask()
#c = RandomDipoles(Rx=40, Ry=40, num=100)
#c.generate_ran()
#d = c.B_z_ran()
#a.convolveFlat(d)

# to convolve a single dipole with a tilted SQUID layout
#example = twoAnglesConvolve()
#example.Convolve()