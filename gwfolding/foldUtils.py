def GPStoGreenwichMeanSiderealTime(GPS, usemod=True):
    #from International Earth Rotation Service, measured in 2002
    #see http://hypertextbook.com/facts/2002/JasonAtkins.shtml
    #wearth=72.92115090*1e-6;
    import numpy as np
    wearth = 2 * np.pi * (1/365.2425 + 1) / 86400
    
    # same thing in hours / sec
    w = wearth / np.pi * 12
    
    # GPS time for 1/1/2000 00:00:00 UTC
    GPS2000=630720013
    # sidereal time at GPS2000,in hours
    # from http://www.csgnetwork.com/siderealjuliantimecalc.html
    sT0=6+39/60+51.251406103947375/3600

    siderealTime=w*(GPS-GPS2000)+sT0
    if usemod:
        siderealTime=np.mod(siderealTime,24)   
    return siderealTime

def exists(val):
    if val in locals() or val in globals():
        return True
    else:
        return False

def smooth(a,WSZ):
    # a: NumPy 1-D array containing the data to be smoothed
    # WSZ: smoothing window size needs, which must be odd number,
    # as in the original MATLAB implementation
    import numpy as np
    if WSZ % 2 == 0:
        WSZ = WSZ - 1
    a = np.reshape(a,a.size)
    out0 = np.convolve(a,np.ones(WSZ,dtype=int),'valid')/WSZ    
    r = np.arange(1,WSZ-1,2)
    start = np.cumsum(a[:WSZ-1])[::2]/r
    stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
    return np.concatenate((  start , out0, stop  ))