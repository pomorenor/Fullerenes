import numpy as np
from math import pi
import matplotlib.pyplot as plt
from math import sqrt
import imageio
from matplotlib import animation as anim

class Parameters:

    def __init__(self, xmax: float, res: int, dt: float,
     timesteps: int, im_time: bool) -> None:

     self.xmax = xmax
     self.res = res
     self.dt = dt
     self.timesteps = timesteps
     self.im_time = im_time

     self.dx = 2*xmax/res
     self.x = np.arange(-xmax +xmax/res, xmax, self.dx)
     self.dk = pi/xmax
     self.k = np.concatenate((np.arange(0, res / 2), np.arange(-res / 2,0)))*self.dk


class Operators:

    def __init__(self, res: int) -> None:

        self.V = np.empty(res, dtype = complex)
        self.R = np.empty(res, dtype = complex)
        self.K = np.empty(res, dtype = complex)
        self.wfc = np.empty(res, dtype = complex)


def init(par: Parameters, voffset:float, wfcoffset:float) -> Operators:

    opr = Operators(len(par.x))
    opr.V = 0
    #opr.V = 0.5 * (par.x - voffset) ** 2
    opr.wfc = np.exp(-((par.x   ) ** 2) / 2, dtype=complex)

    if par.im_time:
        opr.K = np.exp(-0.5*(par.k**2)*par.dt)
        opr.R = np.exp(-0.5*opr.V*par.dt)
    else:
        opr.K = np.exp(-0.5*(par.k**2)*par.dt*1j)
        opr.R = np.exp(-0.5*opr.V*par.dt*1j)
    return opr






def split_op(par: Parameters, opr: Operators) -> None:

    #file_names = []

    # we need not define x
    density_each_time_step = []

    for i in range(par.timesteps):

        opr.wfc *= opr.R

        opr.wfc = np.fft.fft(opr.wfc)

        opr.wfc *= opr.K

        opr.wfc *= opr.R

        density = np.abs(opr.wfc)**2

        density_each_time_step.append(density)


        if par.im_time:
            renorm_factor = sum(density)*par.dx
            opr.wfc /= sqrt(renorm_factor)





    fig, axs = plt.subplots()

    def animate(frame):

        axs.clear()
        axs.set_xlabel("Position")
        axs.set_ylabel("Probabilit Density")
        graph = plt.plot(par.x, density_each_time_step[frame])

        return density_each_time_step


    movie = anim.FuncAnimation(fig, animate, frames = len(range(par.timesteps)), interval = 50 ,repeat = False)
    plt.show()








def main() -> None:
    par = Parameters(5.0, 500, 0.05, 20000, False)

    # Starting wavefunction slightly offset so we can see it change
    opr = init(par, 0.0, -1.00)
    split_op(par, opr)

if __name__ == "__main__":
    main()
