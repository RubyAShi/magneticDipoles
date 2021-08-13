from singleDipole import SingleDipole
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

class RandomDipoles():
    def __init__(self, Rx=10, Ry=10, delta=.1, N=2 * np.float32(1e7),
                 z0=2, phi=0, sigma=0, num=1, yes=1,
                 xString=None, yString=None, angString=None, magni=None):
        self.Rx = Rx
        self.Ry = Ry
        self.delta = delta
        self.N = N
        self.z0 = z0
        self.phi = phi
        self.sigma = sigma
        self.num = num
        self.yes = yes
        self.xString = xString
        self.yString = yString
        self.angString = angString
        self.magni = magni


    # takes in canvas size and number of random dipoles
    # to produce two vectors that contains info of dipole location
    # and one vector that contains dipole tilt angle
    # randomize magnitude of each dipole if yes is 1
    def generate_ran(self):

        Rx = self.Rx
        Ry = self.Ry
        delta = self.delta
        n = self.num
        N = self.N
        yes = self.yes


        del_str = str(delta)
        deci = del_str[::-1].find('.')

        x_cen = np.round(np.random.uniform(-Rx, Rx, size=n), deci)
        y_cen = np.round(np.random.uniform(-Ry, Ry, size=n), deci)
        ang = np.round(np.random.uniform(0, 360, size=n), 3)
        self.xString = x_cen
        self.yString = y_cen
        self.angString = ang
        if yes == 1:
            magni = np.round(np.random.uniform(N / 10, 10 * N, size=n))
        else:
            magni = N * np.ones(n)
        self.magni = magni

    # compute z component of n random dipoles
    # each dipole has an arbitrary center, magnitude and tilt angle

    def B_z_ran(self):
        n = self.num
        Rx = self.Rx
        Ry = self.Ry
        delta = self.delta
        x_cen = self.xString
        y_cen = self.yString
        magni = self.magni
        ang = self.angString
        z = self.z0

        mesh_size = int(2 * Rx / delta) + 1
        B_total = np.zeros([mesh_size, mesh_size])

        for i in range(n):
            single = SingleDipole(Rx, Ry, delta, magni[i], x_cen[i], y_cen[i], z, ang[i])
            B_total += single.Bz_no_plot()
        fig, ax = plt.subplots()
        im = ax.imshow(B_total, cmap=cm.RdYlGn,
                       origin='lower', extent=[-Rx, Rx, -Ry, Ry])
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label('G', size=14)
        plt.xlabel('um', size=14)
        plt.ylabel('um', size=14)
        plt.title('field strength')
        plt.show()

        return B_total

    def B_tilt_two_ran(self):
        n = self.num
        Rx = self.Rx
        Ry = self.Ry
        delta = self.delta
        x_cen = self.xString
        y_cen = self.yString
        magni = self.magni
        ang = self.angString
        z = self.z0

        mesh_size = int(2 * Rx / delta) + 1
        B_total = np.zeros([mesh_size, mesh_size])

        for i in range(n):
            single = SingleDipole(Rx, Ry, delta, magni[i], x_cen[i], y_cen[i], z, ang[i])
            B_total += single.B_tilt_two_no_plot()
        fig, ax = plt.subplots()
        im = ax.imshow(B_total, cmap=cm.RdYlGn,
                       origin='lower', extent=[-Rx, Rx, -Ry, Ry])
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label('G', size=14)
        plt.xlabel('um', size=14)
        plt.ylabel('um', size=14)
        plt.title('field strength')
        plt.show()

        return B_total



