from dipoles import SquidLayout
from dipoles import FlatConvolve
from dipoles import RandomDipoles
from dipoles import SingleDipole
from dipoles import twoAnglesConvolve
import numpy as np

# to generate 3um SQUID layout
#a = SquidLayout()
#b = a.IBM3um_mask()

# to scale a 3um SQUID layout according to aligning angle
#a = SquidLayout()
#b = a.IBM3um_mask()
#b = a.y_scales()

# to generate magnetic field of a number of randomly aligned dipoles
# to change height use z0
#c = RandomDipoles(Rx=40, Ry=40, num=100)
#c.generate_ran()
#d = c.B_z_ran()

# to convolve magnetic field from above with SQUID pickup loop
#c = FlatConvolve(Rx=40, Ry=40, num=1000, N = 2 * np.float32(1e4))
#c.IBM3um_mask()
#c.generate_ran()
#c.ConvolveFlat()


# to convolve a single dipole with a tilted SQUID layout
#example = twoAnglesConvolve()
#example.Convolve()