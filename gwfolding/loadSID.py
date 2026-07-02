def loadSID(SIDFile, GPSStart, segDuration, SigmaCut, badGPS, ifo1, ifo2):
#    
# Load SID from a frame 
#
# [ statistic, invVariance, isStationary ] = ...
#   loadSID (SIDFile, GPSStart, segDuration, SigmaCut, ifo1, ifo2);
#
# SIDFile      : String: File containing the SID
# GPSStart     : Interer: Start time
# segDuration  : Real/Integer: Suration of each segment (possibly 52sec)
# SigmaCut     : Structure: freq bin is tagged non-stationary, if its difference
#                with the local PSD and adjacent PSD differ by more (less) than
#                SigmaCut.maxD (-SigmaCut.minD). SigmaCut.smoothSpan is the
#                number of f-bin to integrate over for Delta Sigma comparison.
#                If it is larger than the total number of f-bin, only the sum
#                of PSD over all frequencies is compared. If it is negative,
#                the so called `sigma` of radiometer point esimate is compared.
#                SigmaCut.maxA is an absolute (non-relative) cut-off on the PSD.
#                SigmaCut.logfp is logfile pointer.
#                (It may be possible to introduce a frequency dependant cut)
# badGPS       : Integer vector: List of bad GPS segments
# ifo1, ifo2   : Strings: Detectors in the baseline
#
# statistic    : Complex vector: CSD
# invVariance  : Real vector: inverse variance of the statistic
#                To be multiplied with winFactor
# isStationary : Integer vector: 1 if freq bin is stationary, 0 otherwise
# misc         : Structure. Return miscellaneous auxiliary data
#
# 
# Author: Sanjit Mitra <sanjit.mitra@ligo.org>
# Translated to Python by Erik Floden <erik.floden@ligo.org>
    import Fr
    import numpy as np
    from . import foldUtils
    misc = {}
    misc['readerror'] = 1
    
    #======================= Load frequency series data ========================
    
    # Main SID elements
    P1 = Fr.frgetvect(SIDFile, ifo1+ ':AdjacentPSD', GPSStart, segDuration)
    P2 = Fr.frgetvect(SIDFile, ifo2+ ':AdjacentPSD', GPSStart, segDuration)
    #CSD = Fr.frgetvect(SIDFile, ifo1+ifo2+ ':CSD', GPSStart, segDuration)[0]
    reCSD = Fr.frgetvect(SIDFile, ifo1+ifo2+ ':ReCSD', GPSStart, segDuration)
    imCSD = Fr.frgetvect(SIDFile, ifo1+ifo2+ ':ImCSD', GPSStart, segDuration)
    CSD = [np.array(np.vectorize(complex)(reCSD[0], imCSD[0]))]
    del reCSD
    del imCSD
    
    nFreqBin = len(CSD[0])
    # Return values
    statistic = CSD
    invVariance = 1.0 / (P1[0] * P2[0]) # So far we did not encounter zeros in P1, P2
    #============================ Data Quality Cuts ============================
    # The cut is not really applied here, only what should be cut is returned
    
    # Absolute sigma cut
    if SigmaCut['maxA'] > 0.0:
        
        # If smooth range is -ve, use standard DSigma cut using integrated 1/(P_1 P_2)
        if SigmaCut['smoothSpan'] < 0:
            
            misc['absSigma'] = np.sqrt(1.0/np.sum(SigmaCut['w8']**2 * invVariance)) #NEEDS TESTING
            
            isQuiet = np.dot(np.ones((nFreqBin,1)) , np.where((misc['absSigma'] <= SigmaCut['maxA']),1,0)) #NEEDS TESTING
            
        # If smoothing window is larger than full f range, disable frequency wise cut
        # These are handwaving formulae, needs to be validated for production runs
        elif SigmaCut['smoothSpan'] > nFreqBin:
            
            misc['absSigma1'] = np.sqrt(1.0/np.sum(SigmaCut['w8']/P1[0]))
            misc['absSigma2'] = np.sqrt(1.0/np.sum(SigmaCut['w8']/P2[0]))

            isQuiet = np.dot(np.ones((nFreqBin,1)) , np.where(\
                (misc['absSigma1'] <= SigmaCut['maxA']) & (misc['absSigma2'] <= SigmaCut['maxA']),1,0))
            
        # If smoothing window is larger than full f range, disable frequency-wise cut
        else:
            isQuiet = np.where((P1[0] <= SigmaCut['maxA']) & (P2[0] <= SigmaCut['maxA']), 1, 0)
    else:
        isQuiet = np.ones((nFreqBin,1))

    if (SigmaCut['maxD'] > 0.0) or (SigmaCut['minD'] > 0.0):
        naiveP1 = Fr.frgetvect(SIDFile,ifo1+':LocalPSD',GPSStart, segDuration)
        misc['naiveP1'] = naiveP1[0]
        naiveP2 = Fr.frgetvect(SIDFile,ifo2+':LocalPSD',GPSStart, segDuration)
        misc['naiveP2'] = naiveP2[0]
        # If smooth range is -ve, use standard DSigma cut using integrated 1/(P_1 P_2)
        if SigmaCut['smoothSpan'] < 0:
            misc['thrSigma'] = np.sqrt(1.0/np.sum(SigmaCut['w8']**2 * invVariance))
            misc['naiSigma'] = np.sqrt(1.0/np.sum(SigmaCut['w8']**2 / (misc['naiveP1']*misc['naiveP2'])))
            ratio =  misc['naiSigma'] / misc['thrSigma']

            if SigmaCut['minD'] > 0.0:
                isDSigmaSmallMin = np.where((misc['naiSigma'] >= np.dot(SigmaCut['minD'] , misc['thrSigma'])),1,0)
            else:
                isDisgmaSmallMin = 1
                
            if SigmaCut['maxD'] > 0.0:
                isDSigmaSmallMax = np.where((misc['naiSigma'] <= np.dot(SigmaCut['maxD'] , misc['thrSigma'])),1,0)
            else:
                isDSigmaSmallMax = 1
            isStationary = np.dot(np.where((isDSigmaSmallMin) & (isDSigmaSmallMax),1,0) , isQuiet)
        # If smoothing window is larger than full f range, disable frequency-wise cut
        elif SigmaCut['smoothSpan'] > nFreqBin:
            print(SigmaCut['smoothSpan'], 'line 110')
            misc['thrSigma1'] = np.sqrt(1.0/np.sum(SigmaCut['w8']/P1[0]))
            misc['thrSigma2'] = np.sqrt(1.0/np.sum(SigmaCut['w8']/P2[0]))
            misc['naiSigma1'] = np.sqrt(1.0/np.sum(SigmaCut['w8']/naiveP1[0]))
            misc['naiSigma2'] = np.sqrt(1.0/np.sum(SigmaCut['w8']/naiveP2[0]))
            
            if SigmaCut['minD'] > 0.0:
                isDSigmaSmallMin = np.where((misc['naiveP1'] >= np.dot(SigmaCut['minD'] , misc['thrSigma1'])) & \
                                            (misc['naiveP2'] >= np.dot(SigmaCut['minD'] , misc['thrSigma2'])), 1, 0)
            else:
                isDSigmaSmallMax = 1
            
            if SigmaCut['maxD'] > 0.0:
                isDSigmaSmallMax = np.where((misc['naiveP1'] <= np.dot(SigmaCut['maxD'] , misc['thrSigma1'])) & \
                                             (misc['naiveP2'] <= np.dot(SigmaCut['maxD'] , misc['thrSigma2'])), 1, 0)
                
            isStationary = np.dot(np.where((isDSigmaSmallMin) and (isDSigmaSmallMax),1,0) , isQuiet)
        # Otherwise use f wise sigma cuts, but smooth over a given range
        # Incorporating weights here is possible, but should be done with care
        else:
            # Compare naive and adjacent PSDs
            if SigmaCut['maxD'] > 0.0:
                #print(naiveP1[0]-P1[0])
                #print(np.dot(SigmaCut['maxD'],P1[0]))
                isDSigmaSmallMax = np.where(((naiveP1[0]-P1[0]) <= np.dot(SigmaCut['maxD'],P1[0])) & \
                                            ((naiveP2[0]-P2[0]) <= np.dot(SigmaCut['maxD'],P2[0])), 1, 0)
            else:
                isDSigmaSmallMax = np.ones((nFreqBin,1))
              
            if SigmaCut['minD'] > 0.0:
                isDSigmaSmallMin = np.where(((P1[0]-naiveP1[0]) >= np.dot(-SigmaCut['minD'],P1[0])) & \
                                            ((naiveP2[0]-P2[0]) >= np.dot(-SigmaCut['minD'],P2[0])), 1, 0)
            else:
                isDSigmaSmallMin = np.ones((nFreqBineqBin,1))
            # Smooth and round the checked values, ignore the end parts
            pyround = np.vectorize(round)
            print(SigmaCut['smoothSpan'], 'line 146')
            isStationary = np.where((pyround(foldUtils.smooth(isQuiet,SigmaCut['smoothSpan']))) & \
                (pyround(foldUtils.smooth(isDSigmaSmallMax,SigmaCut['smoothSpan']))) & \
                (pyround(foldUtils.smooth(isDSigmaSmallMin,SigmaCut['smoothSpan']))),1,0)
            isStationary[0:SigmaCut['smoothSpan']] = 1
            isStationary[SigmaCut['smoothSpan']:] = 1
    else:
 
        # Default is all stationary
        isStationary = isQuiet
    # If bad GPS list was loaded and this segment is one of them
    # (Since number of bad GPS is low, it is still fine to load everything above)
    if GPSStart in badGPS:
        isStationary[:] = 0
        print('WARNING: segment at GPS '+str(GPSStart)+'s was in badGPS list, excluded.',file=SigmaCut['logfp']) 
        
    
    if any(isStationary == 0) and any(isStationary == 1):
        # To reject the whole segment if even one frequency bin is non-stationary, use isStationary[:] = 0
        print('WARNING: Sigma cut applied to some f-bins of the segment at GPS '+str(GPSStart)+'s.',file=SigmaCut['logfp'])
    #================================= Return ==================================
    
    del P1, P2, CSD, isQuiet
    
    if ((SigmaCut['maxD'] > 0.0) or (SigmaCut['minD'] > 0.0)):
        del naiveP1, naiveP2, isDSigmaSmallMax, isDSigmaSmallMin
    
    misc['readerror'] = 0
    
    return statistic, invVariance, isStationary, misc
