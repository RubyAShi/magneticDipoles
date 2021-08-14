import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from PIL import Image, ImageDraw
from scipy.ndimage import convolve

pi = np.pi
# magnetic permearbility in H/m
mu_0 = 4 * pi * np.float32(1e-7)
# Bohr magneton in J/T
m_0 = 9.274 * np.float32(1e-24)

# genetrates magnetic field of a single dipole
class SingleDipole():

    def __init__(self, Rx=10, Ry=10, delta=.1, N=2 * np.float32(1e7), x0=0, y0=0, z0=2, theta=0, phi = 0, sigma = 0):
        # Rx x scanning range in um
        # Ry y scanning range in um
        # delta pixel size
        # z0 scanning height
        # x0, y0 initial dipole location
        # theta dipole angle to x axis in degrees
        # N number of Bohr magnetons
        # phi is aligning angle, 2 to 4 degrees
        # sigma is horizontal tilt or imperfect aligning
        self.Rx = Rx
        self.Ry = Ry
        self.delta = delta
        self.N = N
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.theta = theta
        self.phi = phi
        self.sigma = sigma

    # calculates the z component of the B-filed of a dipole in the xy-plane
    def Bz(self):

        theta_rad = self.theta / 360 * 2 * np.pi
        Rx = self.Rx
        Ry = self.Ry
        delta = self.delta
        x0 = self.x0
        y0 = self.y0
        z0 = self.z0
        N = self.N

        # mesh the plotting area
        x = np.arange(-Rx, Rx + delta, delta)
        y = np.arange(-Ry, Ry + delta, delta)
        X, Y = np.meshgrid(x, y)
        Z = z0 * np.ones((len(x), len(y)))
        X0 = x0 * np.ones((len(x), len(y)))
        Y0 = y0 * np.ones((len(x), len(y)))

        # a position in space. r-vector
        R0 = np.sqrt((X - X0) ** 2 + (Y - Y0) ** 2 + Z ** 2)

        # calculates the magnetic field in Gauses
        m_dot_r = (X - X0) * np.cos(theta_rad) + (Y - Y0) * np.sin(theta_rad)
        B_z = mu_0 / (4 * pi) * N * m_0 * 3 * z0 * m_dot_r / R0 ** 5 * np.float32(1e18) * np.float32(1e4)

        fig, ax = plt.subplots()
        im = ax.imshow(B_z, cmap=cm.RdYlGn,
                       origin='lower', extent=[-Rx, Rx, -Ry, Ry])
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label('G', size=14)
        plt.xlabel('um', size=14)
        plt.ylabel('um', size=14)
        plt.title('field strength')
        plt.show()
        return B_z

    # calculates the z component of the B-filed of a dipole in the xy-plane
    # skip plotting
    def Bz_no_plot(self):


        theta_rad = self.theta / 360 * 2 * np.pi
        Rx = self.Rx
        Ry = self.Ry
        delta = self.delta
        x0 = self.x0
        y0 = self.y0
        z0 = self.z0
        N = self.N

        # mesh the plotting area
        x = np.arange(-Rx, Rx + delta, delta)
        y = np.arange(-Ry, Ry + delta, delta)
        X, Y = np.meshgrid(x, y)
        Z = z0 * np.ones((len(x), len(y)))
        X0 = x0 * np.ones((len(x), len(y)))
        Y0 = y0 * np.ones((len(x), len(y)))

        # a position in space. r-vector
        R0 = np.sqrt((X - X0) ** 2 + (Y - Y0) ** 2 + Z ** 2)

        # calculates the magnetic field in Gauses
        m_dot_r = (X - X0) * np.cos(theta_rad) + (Y - Y0) * np.sin(theta_rad)
        B_z = mu_0 / (4 * pi) * N * m_0 * 3 * z0 * m_dot_r / R0 ** 5 * np.float32(1e18) * np.float32(1e4)
        return B_z

    # calculates the B field projected to a plane slightly tilted to the z-plane
    # there are two tilt angels
    # phi is the aligning angle
    # sigma is a smaller twist of aligning angle
    # plane is characterized by normal vector (sin(\sigma), -cos(\sigma)sin(\phi), cos(\phi)cos(\sigma))
    def B_tilt_two(self):

        theta_rad = self.theta / 360 * 2 * np.pi
        Rx = self.Rx
        Ry = self.Ry
        delta = self.delta
        x0 = self.x0
        y0 = self.y0
        z0 = self.z0
        N = self.N

        phi_rad = self.phi / 360 * 2 * pi
        sigma_rad = self.sigma / 360 * 2 * pi

        # mesh the plotting area
        x = np.arange(-Rx, Rx + delta, delta)
        y = np.arange(-Ry, Ry + delta, delta)
        X, Y = np.meshgrid(x, y)
        Z = z0 * np.ones((len(x), len(y)))
        X0 = x0 * np.ones((len(x), len(y)))
        Y0 = y0 * np.ones((len(x), len(y)))

        # a position in space. r-vector
        R0 = np.sqrt((X - X0) ** 2 + (Y - Y0) ** 2 + Z ** 2)

        # calculates the magnetic field in Gauses
        m_dot_r = (X - X0) * np.cos(theta_rad) + (Y - Y0) * np.sin(theta_rad)

        first = Z * np.cos(phi_rad) * np.cos(sigma_rad)
        - (Y - Y0) * np.sin(phi_rad) * np.cos(sigma_rad) + (X - X0) * np.sin(sigma_rad)

        second = np.sin(theta_rad) * np.sin(phi_rad) * np.cos(sigma_rad) - np.cos(theta_rad) * np.sin(sigma_rad)

        B_tilt = mu_0 / (4 * pi) * N * m_0 * (3 * first * m_dot_r / R0 ** 5 + second / R0 ** 3) * np.float32(
            1e18) * np.float32(1e4)

        fig, ax = plt.subplots()
        im = ax.imshow(B_tilt, cmap=cm.RdYlGn,
                       origin='lower', extent=[-Rx, Rx, -Ry, Ry])
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label('G', size=14)
        plt.xlabel('um', size=14)
        plt.ylabel('um', size=14)
        plt.title('field strength')
        plt.show()
        return B_tilt

    # calculates the same with B_tilt_two but no plot to save simulation time
    def B_tilt_two_no_plot(self):

        theta_rad = self.theta / 360 * 2 * np.pi
        Rx = self.Rx
        Ry = self.Ry
        delta = self.delta
        x0 = self.x0
        y0 = self.y0
        z0 = self.z0
        N = self.N

        phi_rad = self.phi / 360 * 2 * pi
        sigma_rad = self.sigma / 360 * 2 * pi

        # mesh the plotting area
        x = np.arange(-Rx, Rx + delta, delta)
        y = np.arange(-Ry, Ry + delta, delta)
        X, Y = np.meshgrid(x, y)
        Z = z0 * np.ones((len(x), len(y)))
        X0 = x0 * np.ones((len(x), len(y)))
        Y0 = y0 * np.ones((len(x), len(y)))

        # a position in space. r-vector
        R0 = np.sqrt((X - X0) ** 2 + (Y - Y0) ** 2 + Z ** 2)

        # calculates the magnetic field in Gauses
        m_dot_r = (X - X0) * np.cos(theta_rad) + (Y - Y0) * np.sin(theta_rad)

        first = Z * np.cos(phi_rad) * np.cos(sigma_rad)
        - (Y - Y0) * np.sin(phi_rad) * np.cos(sigma_rad) + (X - X0) * np.sin(sigma_rad)

        second = np.sin(theta_rad) * np.sin(phi_rad) * np.cos(sigma_rad) - np.cos(theta_rad) * np.sin(sigma_rad)

        B_tilt = mu_0 / (4 * pi) * N * m_0 * (3 * first * m_dot_r / R0 ** 5 + second / R0 ** 3) * np.float32(
            1e18) * np.float32(1e4)
        return B_tilt

    # calculates the difference betweeb Bz and B_tilt_two
    def get_dif_z(self):

        Rx = self.Rx
        Ry = self.Ry

        B1 = self.Bz_no_plot()
        B2 = self.B_tilt_two_no_plot()
        dif = B2 - B1
        fig, ax = plt.subplots()
        im = ax.imshow(dif, cmap=cm.RdYlGn,
                       origin='lower', extent=[-Rx, Rx, -Ry, Ry])
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label('G', size=14)
        plt.xlabel('um', size=14)
        plt.ylabel('um', size=14)
        plt.title('field strength')
        plt.show()

#genetrates magnetic field of many randomly aligned dipoles
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

class SquidLayout():

    def __init__(self, type = 'IBM3um', xRange=[0,0], yRange=[0,0],z = 2, delta = .1,
                 xCenter = 35, yCenter = 38, mask = None, scaledMask = None, phi = 4, sigma = 0):
        self.type = type
        self.xRange = xRange
        self.yRange = yRange
        self.z = z
        self.delta = delta
        self.xCenter = xCenter
        self.yCenter = yCenter
        self.mask = mask
        self.scaledMask = scaledMask
        self.phi = phi
        self.sigma = sigma

    def get_rect(self, x, y, width, height, angle):
        rect = np.array([(0, 0), (width, 0), (width, height), (0, height), (0, 0)])
        theta = (np.pi / 180.0) * angle
        R = np.array([[np.cos(theta), -np.sin(theta)],
                      [np.sin(theta), np.cos(theta)]])
        offset = np.array([x, y])
        transformed_rect = np.dot(rect, R) + offset
        return transformed_rect

    # SQUID geometry of a 3um IBM SQUID
    # grid size is 0.1
    # outputs a matrix "mask" of size 69 * 69
    def IBM3um_mask(self):
        if self.type != 'IBM3um':
            print('wrong type')
        else:
            data = np.zeros(69 * 69).reshape(69, 69)

            # Convert the numpy array to an Image object.
            img = Image.fromarray(data)

            # Draw a rotated rectangle on the image.
            draw = ImageDraw.Draw(img)
            rect = self.get_rect(x=32, y=0, width=4, height=7, angle=0)
            draw.polygon([tuple(p) for p in rect], fill=1)
            draw.ellipse((4, 8, 64, 68), fill=1)
            # Convert the Image data to a numpy array.
            mask = np.asarray(img)

            # Display the result using matplotlib.
            # this plots from top to bottom
            # first row of matrix is the top



            fig, ax = plt.subplots(1, 2, num = 3)
            ax[0].imshow(mask, cmap=plt.cm.gray)
            ax[0].set_title('matrix plot')
            ax[1].imshow(mask, cmap=plt.cm.gray, extent=[-3.4, 3.4, -3, 3.8])
            ax[1].set_title('real unit plot')
            plt.xlabel("um", labelpad=-2)
            plt.ylabel("um", labelpad=-5)
            plt.show()

            self.type = 'IBM3um'
            self.xRange = [-3.4, 3.4]
            self.yRange = [-3, 3.8]
            self.delta = .1
            self.xCenter = 35
            self.yCenter = 38
            self.mask = mask


    def y_scales(self):
        # add tilt angle info to mask

        type = self.type
        xRange = self.xRange
        yRange = self.yRange
        delta = self.delta
        xCenter = self.xCenter
        yCenter = self.yCenter
        mask = self.mask
        phi = self.phi
        sigma = self.sigma
        z = self.z

        # phi is aligning angle
        phi_rad = phi / 360 * 2 * pi
        sigma_rad = sigma / 360 * 2 * pi

        # mesh the plotting area
        x = np.arange(xRange[0], xRange[1] + delta, delta)
        y = -np.arange(-yRange[1], -yRange[0] + delta, delta)
        X, Y = np.meshgrid(x, y)

        delta_z = -np.sin(sigma_rad) * np.cos(phi_rad) * X + np.sin(phi_rad) * Y
        Z = z * np.ones((len(mask), len(mask[0, :]))) + delta_z
        new_mask = mask * Z

        fig, ax = plt.subplots(2,2, num=2)
        fig.tight_layout()
        im1= ax[0, 0].imshow(delta_z, cmap=plt.cm.gray)
        ax[0, 0].set_title('matrix plot')
        cbar1 = fig.colorbar(im1, ax=ax[0,0])
        cbar1.ax.set_title('um')

        im2 = ax[0, 1].imshow(Z, cmap=plt.cm.gray)
        ax[0, 1].set_title('z height')
        cbar2 = fig.colorbar(im2, ax=ax[0,1])
        cbar2.ax.set_title('um')

        im3 = ax[1, 0].imshow(new_mask, cmap=plt.cm.gray)
        ax[1, 0].set_title('weighted mask')
        cbar3 =fig.colorbar(im3, ax=ax[1, 0])
        cbar3.ax.set_title('um')

        im4 = ax[1, 1].imshow(new_mask, cmap=plt.cm.gray, extent=[-3.4, 3.4, -3, 3.8])
        ax[1, 1].set_title('weighted mask')
        plt.xlabel('um')
        plt.ylabel('um')
        cbar4 = fig.colorbar(im4, ax=ax[1, 1])
        cbar4.ax.set_title('um')
        plt.show()

        self.scaledMask = new_mask


class twoAnglesConvolve(SingleDipole, SquidLayout):
    def __init__(self, Rx=10, Ry=10, B_delta= None, N=2 * np.float32(1e7) , x0=0, y0=0, z0=2, theta=0, phi=4,
                 sigma=2, type='IBM3um', xRange=[0,0], yRange=[0,0], z = 1, delta = .1,
                 xCenter = 5,yCenter = 6,mask = None,scaledMask =None):
        #Rx = 10, Ry = 10, delta = .1, N = 2 * np.float32(1e7), x0 = 0, y0 = 0, z0 = 2, theta = 0, phi = 0, sigma = 0
        #type = 'IBM3um', xRange = [0, 0], yRange = [0, 0], z = 2, delta = .1,
        #xCenter = 35, yCenter = 38, mask = None, scaledMask = None, phi = 4, sigma = 0
        SingleDipole.__init__(self, Rx,Ry,B_delta,N,x0,y0,z0,theta,phi,sigma)
        SquidLayout.__init__(self, type, xRange, yRange, z, delta,
                 xCenter, yCenter, mask, scaledMask)
        self.z = self.z0
        self.B_delta = self.delta


    def Convolve(self):
        #print(self.z0)
        if self.type == 'IBM3um':
            self.IBM3um_mask()
        self.z = self.z0
        #print(self.z)
        self.y_scales()
        x_min = self.xRange[0]
        x_max = self.xRange[1]
        y_min = self.yRange[0]
        y_max = self.yRange[1]
        delta = self.delta
        X = np.arange(x_min, x_max + delta, delta)
        Y = -np.arange(-y_max, -y_min + delta, delta)
        Rx = self.Rx
        Ry = self.Ry
        scanx = np.arange(-Rx, Rx + delta, delta)
        scany = np.arange(-Ry, Ry + delta, delta)
        total_B = np.zeros([len(scanx), len(scany)])
        layout = self.scaledMask
        xcenter = self.xCenter
        ycenter = self.yCenter
        #check_z = np.zeros([len(layout), len(layout[0])])
        #print(self.z0)
        #print(self.delta)
        for x in X:
            for y in Y:
                x0 = round(x, 2)
                y0 = round(y, 2)
                x_position = int(x0 / delta + xcenter - 1)
                y_position = int(-y0 / delta + ycenter)
                # print((x_position, y_position))
                z0 = layout[y_position, x_position]
                # print((x0, y0, x_position, y_position, z0))
                # check_z[x_position, y_position] = layout[x_position, y_position]
                if z0 != 0:
                    self.z0 = z0
                    self.delta = delta
                    self.x0 = -x0
                    self.y0 = -y0
                    total_B += self.B_tilt_two_no_plot()

        Phi_0 = 20.7
        flux = total_B * delta ** 2 / Phi_0 * 1000

        fig, ax = plt.subplots()
        # cmap=cm.YlOrBr
        # cmap=cm.bone
        # cmap=cm.YlGnBu
        cmap = plt.cm.RdYlGn
        im = ax.imshow(flux, cmap, origin='lower', extent=[-Rx, Rx, -Ry, Ry])
        plt.title("magnetic flux simulation")
        plt.xlabel("um")
        plt.ylabel("um")
        cRange = np.max(flux)
        # cRange = 5
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label('m$\Phi_0$', size=14)
        cbar.mappable.set_clim(-cRange, cRange)
        plt.show()

        return flux
        #plt.scatter(scanx, flux[int(Rx / delta) + 1, :], s=5, marker='o')
        #plt.xlabel('um', size=14)
        #plt.ylabel('m$\Phi_0$', size=14)
        #plt.title('cross-section', size=14)
        #plt.show()

        #fig, ax = plt.subplots()
        #im = plt.imshow(check_z, cmap=plt.cm.gray)
        #plt.title("matrix plot")
        #cbar = fig.colorbar(im, ax=ax)
        #cbar.set_label('weighted mask', size=14)
        #plt.show()

        #fig, ax = plt.subplots()
        #im = plt.imshow(check_z, cmap=plt.cm.gray, extent=[-3.4, 3.4, -3.4, 3.4])
        #plt.title("real unit plot")
        #plt.xlabel("um")
        #plt.ylabel("um")
        #cbar = fig.colorbar(im, ax=ax)
        #cbar.set_label('weighted mask', size=14)
        #plt.show()

class FlatConvolve(RandomDipoles, SquidLayout):

    def __init__(self, Rx=10, Ry=10, B_delta= None, N=2 * np.float32(1e7) , z0=2, phi = 0,
                 sigma=2, num=1, yes=1, type='IBM3um', xRange=[0,0], yRange=[0,0], z=None, delta=.1,
                 xCenter = 5, yCenter = 6, mask=None, scaledMask =None):
        #Rx = 10, Ry = 10, delta = .1, N = 2 * np.float32(1e7),
        #z0 = 2, phi = 0, sigma = 0, num = 1, yes = 1,
        #xString = None, yString = None, angString = None, magni = None):
        #type = 'IBM3um', xRange = [0, 0], yRange = [0, 0], z = 2, delta = .1,
        #xCenter = 35, yCenter = 38, mask = None, scaledMask = None, phi = 4, sigma = 0
        RandomDipoles.__init__(self, Rx,Ry,B_delta,N,z0,phi,sigma,num,yes)
        SquidLayout.__init__(self, type, xRange, yRange, z, delta,
                 xCenter, yCenter, mask, scaledMask)
        self.z = self.z0
        self.B_delta = self.delta

    def ConvolveFlat(self):
        type = self.type
        Rx = self.Rx
        Ry = self.Ry
        xRange = self.xRange
        yRange = self.yRange
        delta = self.delta
        xCenter = self.xCenter
        yCenter = self.yCenter
        mask = self.mask
        phi = self.phi
        sigma = self.sigma
        z = self.z

        B_ran = self.B_z_ran()
        flux = convolve(B_ran, mask)
        Phi_0 = 20.7
        flux = flux * delta ** 2 / Phi_0 * 1000

        fig, ax = plt.subplots()
        im = ax.imshow(flux, cmap=plt.cm.YlGnBu,
                       origin='lower', extent=[-Rx,Rx,-Ry,Ry])

        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label('m$\Phi_0$', size=14)
        plt.xlabel('um', size=14)
        plt.ylabel('um', size=14)
        plt.title('magnetic flux')
        plt.show()

        return flux
