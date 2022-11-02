import numpy as np


class Parameters:

    def __init__(self, xmax: float, res: int, dt: float,
     timesteps: int, im_time: bool) -> None:

     self.xmax = xmax
     self.res = res
     self.dt = dt
     self.timesteps = timesteps
     self.im_time = im_time

     self.dx = 2*xmax/res
     self.x = np.arange(0, xmax, self.dx)
     self.dk = pi/xmax
     self.k = np.concatenate((np.arange(0, res / 2), np.arange(-res / 2,0)))*self.dk


class Operators:

    def __init__(self, res: int) -> None:

        self.V = np.empty(res, dtype = complex)
        self.R = np.empty(res, dtype = complex)
        self.K = np.empty(res, dtype = complex)
        self.wfc = np.empty(res, dtype = complex)

def __init__(par: Parameters, voffset:float, wfcoffset:float) -> Operators:

    opr = Operators(len(par.x))
    opr.V = -1/par.x - (1-np.exp(-par.x))**2
    
