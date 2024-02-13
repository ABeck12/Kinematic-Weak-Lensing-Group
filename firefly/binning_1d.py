import numpy as np

def binning_1d(arr, binfac):
    """ Bins up two or three column spectral data by a specified factor. """
    binfac = int(binfac)
    nbins  = len(arr) // binfac
    binned = np.zeros((nbins, 1))
    for i in range(len(binned)):
        spec_slice = arr[i*binfac:(i+1)*binfac]
        binned[i] = np.mean(spec_slice)
    return binned[:,0]