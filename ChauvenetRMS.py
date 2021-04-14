import numpy as np
from scipy.special import erfc


def chauvenet(array):
    """Short summary.
    Determines if each element of the ndarray meets the chauvenet criterion.

    Parameters
    ----------
    array : ndarray
        Data to calculate the chauvenet criterion from

    Returns
    -------
    ndarray
        Boolean array indicating which elements meet the criterion

    """
    mean = np.nanmean(array)          # Mean of incoming array
    stdv = np.nanstd(array)            # Standard deviation
    N = np.size(array)        # Number of non-nan elements in the array
    criterion = 1.0/(2*N)         # Chauvenet's criterion
    d = abs(array-mean)/stdv      # Distance of a value to mean in stdv's
    prob = erfc(d)                # Area normal dist.
    return prob > criterion       # Use boolean array outside this function


def calculatenoise(fullimage):
    """Short summary.
    Determines the noise as the rms of the pixels that satisfy the chauvenet
    criterion. This is to remove the true emission from the background noise
    so as to calculate the rms of the sky only.

    Parameters
    ----------
    fullimage : ndarray
        Array of data, which can be in any dimensions.

    Returns
    -------
    noise : float
        Background noise of the data.
    skyimage : ndarray
        Background emission with the emission detected by the chauvenet
        algorithm masked.

    """
    skyimage = fullimage.copy()
    satisfiescriteria = chauvenet(skyimage)
    # As NaNs do not satistfy the criteria, we force it
    satisfiescriteria[np.where(np.isnan(skyimage))] = True
    while not np.all(satisfiescriteria):
        skyimage[np.where(satisfiescriteria == False)] = np.nan
        satisfiescriteria = chauvenet(skyimage)
        satisfiescriteria[np.where(np.isnan(skyimage))] = True

    noise = np.sqrt(np.nanmean((skyimage)**2))
    return noise, skyimage
