def foldSID(paramsFile, writeFrames, firstSeg=-1, lastSeg=-1):
#
# Main routine to fold Stochastic Intermediate Data (SID) frames to one
# sidereal day. The output can be written as both .mat and/or .gwf formats.
#
# Mathematical details can be found in the LIGO technical document: T0900093
#
#
# [foldedInvCov, foldedStatistic, params, foldParams] = foldSID (paramFile, firstSeg, lastSeg, matFile, writeFrames)
#
# paramsFile      : String. The parameter file
# ------------ Optional parameters for parallelization ------------
# firstSeg        : Integer. First segment # to fold from list
#                   Default -1, equivalent to first available
# lastSeg         : Integer. Last segment # to fold from list
#                   Default -1, equivalent to last available
# writeFrames     : Logical. If exists, frame(s) will be written                   
# -----------------------------------------------------------------
#
# foldedInvCov    : Real matrix (t,f). Inverse covariance of folded statistic
# foldedStatistic : Complex matrix (t,f). Folded statistic
# params          : Structure. Parameters read from the first frame
# foldParams      : Structure. Parameters used for folding
#

#
# Author: Sanjit Mitra <sanjit.mitra@ligo.org>
# Translated to Python by Erik Floden <erik.floden@ligo.org>
#
#============================= Import Packages =============================#
    import warnings
    import sys
    import numpy as np
    import Fr
    import math
    from .loadSID import loadSID
    from .loadSIDParams import loadSIDParams
    from .foldUtils import GPStoGreenwichMeanSiderealTime, exists, smooth
    from .writeFSID import writeFSID
    from scipy.io import savemat
#============================= Read Parameters =============================#
    namesVals = {}
    with open(paramsFile) as f:
        for line in f:
            name, var = line.partition(" ")[::2]
            print(name, var)
            if name.strip() != '%' and name !='#' and var != ' ' and len(name.strip())>0:
                namesVals[name.strip()] = str(var)[:-1]
        f.close()
    SigmaCut ={}
    foldParams = {}
    foldParams['paramsFile'] = paramsFile
    

    for ii in namesVals:
        # Each line of this file should have:
        # Column 1: GPSStart, Column 2: path to the corresponding SID frame
        if ii == 'SIDMetaDataFile':
            foldParams[ii] = namesVals[ii]
            SIDMetaDataFile = foldParams[ii]

        # Segment duration in SID frames (not folded SID frames)
        elif ii == 'segmentDuration' or ii == 'segDuration':
            foldParams['segDuration'] = float(namesVals[ii])
            segDuration = foldParams['segDuration']
            if segDuration <= 0:
                raise ValueError('segDuration must be a positive number')

        # Interferometer 1        
        elif ii == 'ifo1':
            foldParams[ii] = namesVals[ii]
            ifo1 = foldParams[ii]

        # Interferometer 2
        elif ii == 'ifo2':
            foldParams[ii] = namesVals[ii]
            ifo2 = foldParams[ii]

        # (OPTIONAL) Apply Sigma cuts, if necessary
        elif ii == 'maxASigma':
            maxASigma = float(namesVals[ii])
            SigmaCut['maxA'] = maxASigma
        elif ii == 'maxDSigma':
            maxDSigma = float(namesVals[ii])
            SigmaCut['maxD'] = maxDSigma
        elif ii == 'minDSigma':
            minDSigma = float(namesVals[ii])
            SigmaCut['minD'] = minDSigma
        elif ii == 'smoothWinSize4SigmaCut':   # #bin for smoothing
            smoothSpan = int(namesVals[ii])
            SigmaCut['smoothSpan'] = smoothSpan
        elif ii =='sigmaCutWeightFile':       # extra weight for broadband cuts
            sigmaCutWeightFile = foldParams[ii]
            SigmaCut['weightFile'] = sigmaCutWeightFile

        # (OPTIONAL) Do not include bad GPS segments derived externally
        elif ii =='badGPSTimesFile':
            foldParams['badGPSTimesFile']=namesVals[ii]
            badGPSTimesFile=foldParams['badGPSTimesFile']

        # (OPTIONAL) To save channels for compatibility with stochastic.m
        elif ii == 'explicit4isotropic': # Obsolete
            foldParams['backwardCompatible'] = True
            backwardCompatible = foldParams['backwardCompatible']
            warnings.warn('Option explicit4isotropic is depreciated, use backwardCompatible')

        # (OPTIONAL) To save channels for compatibility with stochastic.m
        elif ii == 'backwardCompatible':
            foldParams[ii] = int(namesVals[ii])
            backwardCompatible = foldParams[ii]

        # (OPTIONAL) Assume that noise PSDs of the neighbors are very close
        elif ii == 'identicalNeighbors':
            foldParams[ii] = int(namesVals[ii])
            identicalNeighbors = foldParams[ii]

        # (OPTIONAL) Create a jobfile for folded data
        elif ii == 'FSIDJobFile':
            foldParams[ii] = namesVals[ii]
            FSIDJobFile = foldParams[ii]

        # (OPTIONAL) Prefix for FSID frames (suffixes are sidereal time, duration, ..)
        # Frames will be written iff this parameter is supplied
        elif ii == 'FSIDFramePrefix':
            foldParams[ii] = namesVals[ii]
            FSIDFramePrefix = foldParams[ii]

        # (OPTIONAL) Offset the sidereal GPSStart, tofool the pipeline!
        elif ii == 'FSIDGPSOffset':
            FSIDGPSOffset = float(namesVals[ii])
            FSIDGPSOffsetReq = float(namesVals[ii])
            foldParams[ii] = FSIDGPSOffsetReq # May be changed later

        # (OPTIONAL) Correct GPS offset to match the sidereal zero
        elif ii == 'FSIDGPSOffsetCorrect':
            FSIDGPSOffsetCorrect = float(namesVals[ii])
            foldParams[ii] = FSIDGPSOffsetCorrect

        # (OPTIONAL) Maximum number of segments per file, default 30 (26min)
        # If this is set to -1, whole folded data will be written to one frame
        elif ii == 'segmentsPerFrame':
            foldParams[ii] = round(float(namesVals[ii]))
            segmentsPerFrame = foldParams[ii]
            if segmentsPerFrame < 1:
                raise ValueError('segmentsPerFrame must be a positive integer')

        # (OPTIONAL) Version of the Folded SID frame
        elif ii == 'version':
            foldParams[ii] = float(namesVals[ii])
            FSIDVersion = foldParams[ii]

        # (OPTIONAL) Sidereal day in seconds, default 86164
        elif ii == 'siderealDay':
            foldParams[ii] = float(namesVals[ii])
            siderealDay = foldParams[ii]
            if siderealDay <= 0:
                raise ValueError('siderealDay must be a positive number')

        # (OPTIONAL/DEBUG) Overlapping window correction, default true
        elif ii == 'ovlWinCorrection':
            foldParams[ii] = float(namesVals[ii])
            ovlWinCorrection = foldParams[ii]

        # (OPTIONAL) Log file, default stdout
        elif ii == 'logFile':
            foldParams[ii] = namesVals[ii]
            logFile = foldParams[ii]

        # (OPTIONAL) After this many SID frames progress will be written, default 1000
        elif ii == 'logStep':
            foldParams[ii] = float(namesVals[ii])
            logStep = foldParams[ii]
            if logStep < 1:
                raise ValueError('logStep must be a positive integer')

        # (OPTIONAL) Write all, basic output will be written by default
        elif ii == 'verbose':
            foldParams[ii] = float(namesVals[ii])
            verbose = foldParams[ii]

        # Default warning
        else:
            raise ValueError('WARNING: Unknown parameter: ' + ii)

#========================= Check Compulsory Inputs =========================#
    try:
        SIDMetaDataFile
    except NameError:
        raise NameError('SIDMetaDataFile not found in file: '+paramsFile)

    try:
        segDuration
    except NameError:
        raise NameError('segmentDuration not found in file: '+paramsFile)    

#=============================== Log details ===============================#

    try:
        logFile
    except NameError:
        logfp = sys.stdout # stdout
    else:
        logfp = open(logFile, "w")
    foldParams['logfp'] = logfp

    try:
        logStep
    except NameError:
        logStep = 1000 # Write status after processing this many frames
    foldParams['logStep'] = logStep

    try:
        verbose
    except NameError:
        verbose = False
    else:
        verbose = True
#======== Load and sort the list of GPS segments and the file names ========#

# Sorting can be slow. Most likely the list is already sorted, check it first
# If any other data query tool is used, modify this part

    listOfGPSStart = []
    listOfSIDFile = []
    with open(SIDMetaDataFile) as f:
        for line in f:
            GPSStart, SIDFile = line.partition("\t")[::2]
            listOfGPSStart.append(float(GPSStart))
            listOfSIDFile.append(SIDFile[:-1])
        f.close()
    if not (all(listOfGPSStart[i] <= listOfGPSStart[i + 1] for i in range(len(listOfGPSStart)-1))): #if listOfGPSStart is not sorted
        sorted_idx = sorted(range(len(listOfGPSStart)),key=lambda x:listOfGPSStart[x])
        listOfSIDFile = [listOfSIDFile[i] for i in sorted_idx]
        listOfGPSStart = [listOfGPSStart[i] for i in sorted_idx]
        del sorted_idx
        #print('was not sorted, but now is')

    nSIDSeg = len(listOfGPSStart)
    print('\nTotal number of segments available: '+ str(nSIDSeg) ,file=logfp)
#=========================== Default Parameters ============================#
    try:
        maxASigma
    except NameError:
        SigmaCut['maxA'] = -1.0 # i.e. disable Absolute Sigma cut

    try:
        maxDSigma
    except NameError:
        SigmaCut['maxD'] = -1.0 # i.e. disable Delta Sigma cut

    try:
        minDSigma
    except NameError:
        SigmaCut['minD'] = -1.0 # i.e. disable Delta Sigma cut

    try:
        smoothSpan
    except NameError:
        SigmaCut['smoothSpan'] = 1000 # Smooth PSDs before checking stationarity

    if (SigmaCut['maxD'] <= 0.0) and (SigmaCut['minD'] <= 0.0):
        print('\n*** WARNING: NOT applying Delta Sigma cuts ***\n' ,file=logfp)
    else:
        print('\nDiscarding segments outside DSigma ratio ['+str(SigmaCut['minD'])+','+str(SigmaCut['maxD'])+']\n' ,file=logfp)

    if (SigmaCut['maxA'] <= 0.0):
        print('\n*** WARNING: NOT applying absolute Sigma cuts ***\n' ,file=logfp)
    else:
        print('\nDiscarding segments if Sigma exceeds '+str(SigmaCut['maxA'])+'\n' ,file=logfp)

    try:
        sigmaCutWeightFile
    except NameError:
        SigmaCut['weightFile'] = ''

    SigmaCut['logfp'] = logfp

    if (SigmaCut['maxD'] > 0.0) or (SigmaCut['minD'] > 0.0):
        SigmaCut['applyDS'] = True
    else:
        SigmaCut['applyDS'] = False

    if (SigmaCut['maxA'] > 0.0):
        SigmaCut['applyAS'] = True
    else:
        SigmaCut['applyAS'] = False

    sigmaCutData = {}
    sigmaCutData['badSegment'] = []

    try:
        badGPSTimesFile
    except NameError:
        badGPSTimesFile = ''
    else:
        print('\nReading bad GPS list from file '+ badGPSTimesFile+'\n' ,file=logfp)

    try:
        backwardCompatible
    except NameError:
        foldParams['backwardCompatible'] = False
        backwardCompatible = foldParams['backwardCompatible']

    try:
        identicalNeighbors
    except NameError:
        foldParams['identicalNeighbors'] = False
        identicalNeighbors = foldParams['identicalNeighbors']

    if identicalNeighbors:
        print('\n*** ASSUMING neighboring segments have the same noise PSD ***\n' ,file=logfp)

    try:
        ovlWinCorrection
    except NameError:
        ovlWinCorrection = 1
        foldParams['ovlWinCorrection'] = ovlWinCorrection
    if ovlWinCorrection == 0:
        print('\n*** WARNING: NOT correcting for overlapping windows ***\n' ,file=logfp)

    try:
        siderealDay
    except NameError:
        siderealDay = 86164.1 ## SHOULD IT BE A COMPULSORY VARIABLE??
        foldParams['siderealDay'] = siderealDay

    try:
        firstSeg
    except NameError:
        firstSeg = 0
    else:
        if firstSeg == -1:
            firstSeg = 0
        elif type(firstSeg) is str:
            if len(firstSeg) == 0:
                firstSeg = 0
            else: firstSeg = round(float(firstSeg))
    if ((firstSeg < 0) or (firstSeg > nSIDSeg)):
        raise ValueError('first segment is out of range, should be in the range [0,'+str(nSIDSeg)+'-1]')
    foldParams['firstSeg'] = firstSeg

    try:
        lastSeg
    except NameError:
        lastSeg = nSIDSeg-1
    else:
        if lastSeg==-1:
            lastSeg = nSIDSeg-1
        elif type(lastSeg) is str:
            if len(lastSeg) == 0:
                lastSeg = nSIDSeg-1
            else: lastSeg = round(float(lastSeg))
    if ((lastSeg < 0) or (lastSeg > (nSIDSeg-1))):
        raise ValueError('last segment is out of range, should be in the range [0,'+str(nSIDSeg)+'-1]')
    foldParams['lastSeg'] = lastSeg

    if (firstSeg > lastSeg):
        raise ValueError('first segment number should be less than last segment number')

    nSIDSegments = lastSeg - firstSeg + 1

    print('First segment to fold: '+str(firstSeg), file=logfp)
    print('Last segment to fold: '+str(lastSeg), file = logfp)
    print('=> Number of segments to be folded: '+str(nSIDSegments), file=logfp)
    
    try:
        FSIDFramePrefix
    except NameError:
        FSIDFramePrefix = 'DNE'
    
    try: 
        segmentsPerFrame
    except NameError:
        segmentsPerFrame = 'DNE'
        
    if FSIDFramePrefix!='DNE' and segmentsPerFrame == 'DNE':
        segmentsPerFrame = 1 
        foldParams['segmentsPerFrame'] = segmentsPerFrame
    
    try:
        FSIDGPSOffsetCorrect
    except NameError:
        FSIDGPSOffsetCorrect = 1
        foldParams['FSIDGPSOffsetCorrect'] = 1
        
    try:
        FSIDGPSOffset
    except NameError:
        FSIDGPSOffset = 'DNE'
        
    if FSIDFramePrefix != 'DNE':
        if FSIDGPSOffset != 'DNE' and FSIDGPSOffsetCorrect == 1:
            foldParams['FSIDGPSOffset'] = FSIDGPSOffsetReq + siderealDay - \
            siderealDay * GPStoGreenwichMeanSiderealTime(FSIDGPSOffsetReq)/24.0
        else:
            FSIDGPSOffsetReq = 0
            foldParams['FSIDGPSOffset'] = 0
            foldParams['FSIDGPSOffsetCorrect'] = 0
            
    try:
        writeFrames
    except NameError:
        writeFrames = 'DNE'
        
    if FSIDFramePrefix != 'DNE':
        if writeFrames == 'DNE':
            print('WARNING: FSID frames will not be written. However, FSIDFramePrefix',file=logfp)
            print('WARNING: variable will be stored in [mat?] file for combine routine.', file=logfp)
        elif writeFrames:
            print('FSID frames will be written',file=logfp)
            print('GPS offset for FSID frames requested: '+str(FSIDGPSOffsetReq),file=logfp)
            print('Actual GPS offset for FSID frames will be: '+str(foldParams['FSIDGPSOffset'])+'sec', file=logfp)
    else:
        if writeFrames != 'DNE' and writeFrames:
            raise ValueError('FSIDFramePrefix is not set, though writeFrames is selected')

#========================== Free some memory here ==========================#

# SID MetaData can take huge amount of memory, erase the unneccessary ones asap
# Keep some of the boundary elements, just in case

    nExtraMetaDataElem2keep = 1 # MUST BE >=1
    if (firstSeg > nExtraMetaDataElem2keep):
        for i in range(firstSeg - 1 - nExtraMetaDataElem2keep):
            listOfSIDFile[i] = [1]
    if (lastSeg < (nSIDSeg-1) - nExtraMetaDataElem2keep):
        for i in range(lastSeg+nExtraMetaDataElem2keep,len(listOfSIDFile)):
            listOfSIDFile[i] = [1]

#====== Load parameters & data from the 1st SID frame to be analyzed =======#
    SIDFile = listOfSIDFile[firstSeg]
    print(SIDFile)
    params = loadSIDParams(SIDFile, ifo1, ifo2, segDuration)

    # Sanity check
    if (params['segmentDuration'] != segDuration):
        raise ValueError('Segment duration mismatch in first frame '+str(SIDFile));

    # Parameters read from the first SID frame
    fLow = params['flow']
    fHigh = params['fhigh']
    deltaF = params['deltaF']

    ################## FIX FACTORS HERE #########################
    winFactor = params['w1w2bar'] * params['w1w2bar'] / params['w1w2squaredbar']
    winRatio = 0.5 * params['w1w2ovlsquaredbar'] / params['w1w2squaredbar']
    ################## FIX FACTORS HERE #########################

    ### DISABLE OVERLAPPING WINDOW CORRECTION, IF FORCED ###
    if ovlWinCorrection == 0:
        winRatio = 0.0
        winFactor = 1.0
    ### DISABLE OVERLAPPING WINDOW CORRECTION, IF FORCED ###

    # Some extra parameters are inserted in the SID params for easy saving!
    params['winFactor'] = winFactor;
    params['winRatio']  = winRatio;

    print('\nFrom the first SID frame:', file=logfp)
    print('fLow = '+str(fLow)+', fHigh = '+str(fHigh)+', deltaF = '+str(deltaF),file=logfp)
    print('winRatio = '+str(winRatio)+', winFactor = '+str(winFactor), file=logfp)

    # Sigma cut related parameters
    if SigmaCut['applyDS']:
        if SigmaCut['smoothSpan'] < 0:
            print('Applying DSigma cut using integrated weight/(P_1P_2)',file=logfp)
        elif (SigmaCut['smoothSpan'] >= np.floor((fHigh-fLow)/deltaF)):
            print('Applying DSigma cut using integrated PSD', file=logfp)
        else:
            print('Applying f wise DSigma cut by smoothing PSD over '+str(SigmaCut['smoothSpan'])+' bins',file=logfp)

    # Externally provided bad GPS list
    if (badGPSTimesFile) != '':
        badGPS = []
        with open(badGPSTimesFile, "r") as my_file:
            contents = my_file.read()
        temp = contents.split('\n')
        for i in range(len(temp)):
            if temp[i] != '':
                badGPS.append(float(temp[i]))
    else:
        badGPS=[]
    if len(SigmaCut['weightFile'])==0:
        SigmaCut['w8'] = 1
    else:
        sigmaCutW8 = np.fromfile(SigmaCut['weightFile'],dtype=float)
        # Linear interpolation, zero outside interval for easy band limiting of sums
        # Log interpolation messes up zeros
        SigmaCut['w8'] = sp.interpolate.interp1d(sigmaCutW8[:,0],sigmaCutW8[:,1],kind='linear')(np.arange(fLow,fHigh+1,deltaF))
        print('Using weights by interpolating data to apply Sigma Cuts from:\n'+ SigmaCut['weightFile'],file=logfp)
        print('Frequency range considered: ['+str(np.min(sigmaCutW8[:,0]))+','+str(np.max(sigmaCutW8[:,1]))+']',file=logfp)
        del SigmaCutW8
    # Save the final set of data quality cut parameters
    foldParams['SigmaCut'] = SigmaCut
#============================ Output parameters ============================#

    foldedSegDuration = segDuration/2 # Factor of 2 to account for 50% overlap
    nFreqBin = int((fHigh - fLow) / deltaF + 1e-9) + 1
    #nFreqBin = math.floor((fHigh-fLow)/deltaF) + 1 #!= ceil if (fHigh-fLow)%deltaF==0
    nSegment = round(siderealDay/foldedSegDuration)
    foldParams['segDist'] = np.zeros((nSegment,1))

    params['foldedSegmentDuration'] = foldedSegDuration
    foldParams['foldedSegDuration'] = foldedSegDuration
    foldParams['nSegment'] = nSegment
    foldParams['nFreqBin'] = nFreqBin

    print('Number of folded segments will be: '+str(nSegment),file=logfp)
    print('Number of f bin in each segment will be: '+str(nFreqBin),file=logfp)

    print('Estimated volume of folded data (memory & disk): '+str(math.ceil(40*nFreqBin*nSegment/1024/1024))+'MB',file=logfp)
    if backwardCompatible:
        print('Extra diskspace for backward compatibility: '+str(math.ceil(24*nFreqBin*nSegment/1024/1024))+'MB',file=logfp)

    print('Total memory & disk usage will be more depending on other parameters.',file=logfp)

    # Allocate memory for output
    foldedInvCov = {}
    foldedStatistic = np.zeros((nSegment,nFreqBin),dtype=complex)
    foldedInvCov['Diag'] = np.zeros((nSegment,nFreqBin))
    foldedInvCov['Prev'] = np.zeros((nSegment,nFreqBin))
    foldedInvCov['Next'] = np.zeros((nSegment,nFreqBin))

    # Store the so called naive and theoretical sigmas
    if SigmaCut['applyDS'] and (SigmaCut['smoothSpan'] < 0):
        sigmaCutData['GPSStart'] = np.zeros((nSIDSegments,1))
        sigmaCutData['naiSigma'] = np.zeros((nSIDSegments,1))
        sigmaCutData['thrSigma'] = np.zeros((nSIDSegments,1))

#================================ MAIN LOOP ================================


    # Before entering the main loop, read the first frame
    GPSStart = listOfGPSStart[firstSeg]
    SIDFile = listOfSIDFile[firstSeg]
    statistic, invVariance, isStationary, misc = \
        loadSID(SIDFile, GPSStart, segDuration, SigmaCut, badGPS, ifo1, ifo2)

    varSigmaSqInv = winFactor * invVariance

    if not any(isStationary):
        sigmaCutData['badSegment'].append(GPSStart)

    params['GPSFirstStart'] = GPSStart

    # Store the so called the naive and theoretical sigmas
    # Could move inside the main loop to have one such code block, but this is safer
    if SigmaCut['applyDS'] and (SigmaCut['smoothSpan'] < 0):
        sigmaCutData['GPSStart'][0] = GPSStart
        sigmaCutData['naiSigma'][0] = misc['naiSigma']
        sigmaCutData['thrSigma'][0] = misc['thrSigma']

    # Also, load the previous SID frame, if it exists
    if firstSeg > 0:

        GPSStartPrev = listOfGPSStart[firstSeg-1]
        SIDFile = listOfSIDFile[firstSeg-1]
        statisticPrev, invVariance, isStationaryPrev = \
            loadSID(SIDFile, GPSStartPrev, segDuration, SigmaCut, badGPS, ifo1, ifo2)
        varSigmaSqInvPrev = winFactor * invVariance
    else:

        # If the first frame analyzed is the first frame available
        # set everything to zero, maintaining proper dimensions
        GPSStartPrev = 0
        statisticPrev = 0.0 * statistic[0]
        varSigmaSqInvPrev = 0.0 * varSigmaSqInv
        isStationaryPrev = 1 + 0 * isStationary

    # Account for offset between siderealStart of FSID with GPSStart of SID frames
    sumStartOffset = 0.0

    # Loop over the SID segments
    print('Reading and folding SID frames...',file=logfp)

    for iSIDSeg in range(nSIDSegments):

        thisSeg = iSIDSeg + firstSeg

        # Corresponding folded (sidereal day) segment
        # NOTE: floor(x) + 1 != ceil(x) when x is an integer
        siderealTime = GPStoGreenwichMeanSiderealTime(GPSStart)
        actualFoldSeg = siderealTime/24.0*nSegment
        foldSeg = round(actualFoldSeg)
        sumStartOffset = sumStartOffset + (actualFoldSeg-foldSeg)*foldedSegDuration
        foldSeg = np.mod(foldSeg,nSegment)

        # Count number of good segments being folded in a given segment
        if np.any(isStationary):
            foldParams['segDist'][foldSeg] = foldParams['segDist'][foldSeg] + 1

        # Show which segment is going where
        if verbose:
            print('GPS = '+str(GPSStart)+', sidereal time = '+str(siderealTime)+\
                ', folded segment # = '+str(foldSeg),file=logfp)

        # Load the NEXT segment here, unless the current segment is the last
        if thisSeg < (nSIDSeg-1):
            GPSStartNext = listOfGPSStart[thisSeg+1]
            SIDFile = listOfSIDFile[thisSeg+1]
            statisticNext, invVariance, isStationaryNext, miscNext = \
                loadSID(SIDFile, GPSStartNext, segDuration, SigmaCut, \
                badGPS, ifo1, ifo2)
            varSigmaSqInvNext = winFactor * invVariance
            if not np.any(isStationaryNext):
                sigmaCutData['badSegment'].append(GPSStartNext)

            # Store the so-called naive and theoretical sigmas
            if SigmaCut['applyDS'] and (SigmaCut['smoothSpan']<0):
                sigmaCutData['GPSStart'][iSIDSeg+1] = GPSStartNext
                sigmaCutData['naiSigma'][iSIDSeg+1] = miscNext['naiSigma']
                sigmaCutData['thrSigma'][iSIDSeg+1] = miscNext['thrSigma']

        else:
        
            # If the last frame analyzed in the last frame available
            # set everything to zero, maintaining proper dimensions
            GPSStartNext = 0
            statisticNext = 0.0 * statistic[0]
            varSigmaSqInvNext = 0.0 * varSigmaSqInv
            isStationaryNext = 1 + 0 * isStationary

        # If previous segment exists and adjacent (allow 1 sec offset)
        prevExists = winRatio * \
            ((thisSeg>1) and ((GPSStart - GPSStartPrev) < foldedSegDuration+1))
        # If next segment exists and adjacent (allow 1 sec offset)
        nextExists = winRatio * \
            ((thisSeg < (nSIDSeg-1)) and ((GPSStartNext - GPSStart) < foldedSegDuration+1))
        
        #======================= Folding Here : Start =============================%
        if not identicalNeighbors: # the actual scenario
            foldedInvCov['Diag'][foldSeg,:] = \
                foldedInvCov['Diag'][foldSeg,:] + isStationary.T * varSigmaSqInv.conj().T

            foldedInvCov['Prev'][foldSeg,:] = foldedInvCov['Prev'][foldSeg,:] + \
                np.dot((0.5 * prevExists) , isStationary.T) * (varSigmaSqInv + varSigmaSqInvPrev).conj().T

            foldedInvCov['Next'][foldSeg,:] = foldedInvCov['Next'][foldSeg,:] + \
                np.dot((0.5 * nextExists) , isStationary.T) * (varSigmaSqInv + varSigmaSqInvNext).conj().T

            foldedStatistic[foldSeg,:] = foldedStatistic[foldSeg,:] + \
                isStationary.T * (varSigmaSqInv * statistic - (\
                np.dot((0.5*prevExists),(varSigmaSqInv + varSigmaSqInvPrev)) * statisticPrev + \
                np.dot((0.5*nextExists),(varSigmaSqInv + varSigmaSqInvNext)) * statisticNext)).conj()

        else: # To match the current stochastic SpH pipeline
            foldedInvCov['Diag'][foldSeg,:] = foldedInvCov['Diag'][foldSeg,:] + \
                np.dot((1.0 - prevExists - nextExists) , isStationary.T) * varSigmaSqInv.conj().T

            foldedStatistic[foldSeg,:] = foldedStatistic[foldSeg,:] + \
                (np.dot((1.0 - prevExists - nextExists) , isStationary.T) * (varSigmaSqInv * statistic).conj())#.T)#.reshape(-1)
        
       #======================= Folding Here : End ===============================#


        # Now shift data elements to process the next segment
        # Current -> Previous
        GPSStartPrev      = GPSStart
        statisticPrev     = statistic
        varSigmaSqInvPrev = varSigmaSqInv
        isStationaryPrev  = isStationary # Redundant

        # Next -> Current
        GPSStart      = GPSStartNext
        statistic     = statisticNext
        varSigmaSqInv = varSigmaSqInvNext
        isStationary  = isStationaryNext
        if (np.remainder(iSIDSeg,logStep) == 0):
            print(str(iSIDSeg)+' of '+str(nSIDSegments)+' done.',file=logfp)



    params['GPSLastEnd'] = GPSStartPrev + foldedSegDuration

    foldParams['startOffset'] = sumStartOffset / nSIDSegments       

#======================= Write folded data to frames =======================#

    print('Writing folded SID frame files to '+ FSIDFramePrefix,file=logfp)
    writeFSID(foldParams, params, foldedInvCov, foldedStatistic)
    print('DONE',file=logfp)
        
#================================= The End =================================#
    try:
        logFile
    except NameError:
        logFile = 'DNE'
    
    if logFile != 'DNE':
        logfp.close()
        

    print('\n** Success: SID from GPS time '+str(params['GPSFirstStart'])+' to '+str(params['GPSLastEnd'])+' folded to 1 sidereal day **\n')
    return foldedInvCov, foldedStatistic, params, foldParams
