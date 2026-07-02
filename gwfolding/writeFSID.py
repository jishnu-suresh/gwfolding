def writeFSID(foldParams, params, foldedInvCov, foldedStatistic):

# Write Folded Stochastic Intermediate Data to frames (.gwf)


# writeFSID (foldParams, params, foldedInvCov, foldedStatistic);

# foldParams      : Structure. Essential elements:
#                    ifo1              : String. Detector 1
#                    ifo2              : String. Detector 2
#                    foldedSegDuration : Real. Duration of folded segments
#                                        This is the foldedSegDuration
#                    nSegment          : Integer. Number of segments to save
#                                        (sidereal day / foldedSegDuration)
#                    segDist           : Integer Array. Number of SID segments
#                                        folded to each sidereal segment
#                    segmentsPerFrame  : Integer. Max # segments per frame
#                                        If -1, put all FSID in one frame
#                    FSIDFramePrefix   : String. Prefix for frame files

#                   Optional Elements:
#                    FSIDGPSOffset     : Integer. Offset to be added to GPSStart
#                    startOffset       : Real. Offset added to siderealStart
# params          : Structure. Parameters from the first SID frame and some more
# foldedInvCov    : Real matrix (t,f). Inverse variance of folded statistic
# foldedStatistic : Complex matrix (t,f). Folded statistic



# Author: Sanjit Mitra <sanjit.mitra@ligo.org>
#         Translated to Python by Erik Floden <erik.floden@ligo.org>


#========= Extract useful parameters from the parameter structures =========#
    import numpy as np
    import sys
    import Fr
    ifo1              = foldParams['ifo1']
    ifo2              = foldParams['ifo2']
    nSegment          = foldParams['nSegment']
    foldedSegDuration = foldParams['foldedSegDuration']
    FSIDFramePrefix   = foldParams['FSIDFramePrefix']
    segmentsPerFrame  = foldParams['segmentsPerFrame']
    FSIDGPSOffset     = foldParams['FSIDGPSOffset']
    if segmentsPerFrame == -1:
        segmentsPerFrame = nSegment

    try:
        backwardCompatible = foldParams['backwardCompatible']
        winFactor = params['winFactor']
    except:
        backwardCompatible = False
        
    try:
        identicalNeighbors = foldParams['identicalNeighbors']
    except:
        identicalNeighbors = False
        
    try:
        SigmaCut = foldParams['SigmaCut']
    except:
        SigmaCut['maxA'] = -1.0
        SigmaCut['maxD'] = -1.0
        SigmaCut['minD'] = -1.0
        SigmaCut['smoothSpan'] = 0
        SigmaCut['logfp'] = sys.stdout
        SigmaCut['applyAS'] = False
        SigmaCut['applyDS'] = False
        
    try:
        FSIDJobFilePtr = open(foldParams['FSIDJobFile'], 'w')
        writeFSIDJobFile = True
    except:
        writeFSIDJobFile = False

    try:
        FSIDGPSOffset = foldParams['FSIDGPSOffset']
    except:
        FSIDGPSOffset = 0

    try:
        siderealDay = foldParams['siderealDay']
    except:
        siderealDay = 86164.1

    try:
        startOffset = foldParams['startOffset']
    except:
        startOffset = 0

    try:
        totNSegments = np.sum(foldParams['segDist'])
    except:
        totNSegments = 0

    try:
        versionString = 'v'+str(foldParams['version'])
    except:
        versionString = ''

    try:
        verbose = foldParams['verbose']
    except:
        verbose = False

    try:
        logfp = foldParams['logfp']
    except:
        logfp = sys.stdout

    # If overlapping window is on, skip half segment to match mid-segment of SIDs
    if foldParams['ovlWinCorrection']:
        startOffset = startOffset + foldedSegDuration/2.0
        
    # Missing segments
    missingSegments = np.where(foldParams['segDist']!=[0])
    if (len(missingSegments) > 0) and (verbose):
        print(missingSegments)
        
    # Entry for job file
    jobCount = 0
    jobStart = -1
    
    #=============== Distribute and write folded data in frames ================#
    
    for iSegment in range(0,nSegment, segmentsPerFrame):
        
        # Number of segments in current frame
        nSegment_thisFile = np.min([nSegment-iSegment, segmentsPerFrame])
        
        # Find missing segments for this file
        missingSegments = np.where(foldParams['segDist'][iSegment:iSegment+nSegment_thisFile]==0)[0]

        # If all segments are missing, write jobfile entry and don't do anything else
        if len(missingSegments) == nSegment_thisFile:
            
            # Write jobfile entry
            if (writeFSIDJobFile and (jobStart >= 0)):
                jobCount = jobCount + 1
                FSIDJobFilePtr.write(str(jobCount)+' '+str(jobStart)+' '+str(jobEnd)+' '+str(jobEnd - jobStart)+'\n')
                jobStart = -1
            
            if verbose:
                print('No data in sidereal time segment '+str(iSegment)+', skipping', file=logfp)
            
            continue
        
        # Start and end sidereal time of data in current frame
        siderealStart = startOffset + (iSegment)*foldedSegDuration
        siderealEnd = siderealStart + foldedSegDuration*nSegment_thisFile
        
        # Every time there is a break in jobfile, move to the next sidereal day
        pyround = np.vectorize(round)
        thisGPSStart = pyround(siderealStart + FSIDGPSOffset + jobCount*siderealDay)
        thisGPSEnd = pyround(siderealEnd + FSIDGPSOffset + jobCount*siderealDay)
        # To disable moving next contiguous segment to next sidereal day
        #thisGPSStart  = pyround(siderealStart + FSIDGPSOffset)
        #thisGPSEnd    = pyround(siderealEnd   + FSIDGPSOffset)
        
        # Entry for job file
        # Keep enough gap before and after a job interval, including buffer
        if writeFSIDJobFile:
            if (jobStart < 0):
                jobStart = thisGPSStart - pyround(3*foldedSegDuration)
            jobEnd = thisGPSEnd + pyround(3*foldedSegDuration)
        
        # Frame File name
        frameFile = (FSIDFramePrefix+ifo1[0]+ifo2[0]+'-FSID'+versionString+'_'+ifo1+ifo2 \
                          +'_'+str(pyround(1000*params['deltaF']))+'mHz-'+str(thisGPSStart)+'-'+str(thisGPSEnd-thisGPSStart)+'.gwf')
        if (verbose):
            print('Writing segments '+str(iSegment)+' through '+str(iSegment+nSegment_thisFile-1)+' to file '+frameFile,file=logfp)

        #=============================== T-F maps ================================#
        
        # Folded statistic
        write_buffer_cplx = foldedStatistic[iSegment:iSegment+nSegment_thisFile,:].conj().T
        if np.isreal(write_buffer_cplx).any():
            write_buffer_cplx = write_buffer_cplx.astype(complex)#np.complex(write_buffer_cplx)
        
        # Make frames directly compatible to isotropic stochastic search, if needed
        
        if (backwardCompatible):
            A = winFactor
            B = foldedInvCov['Diag'][iSegment:iSegment+nSegment_thisFile,:]
            C = foldedInvCov['Prev'][iSegment:iSegment+nSegment_thisFile,:]
            D = foldedInvCov['Next'][iSegment:iSegment+nSegment_thisFile,:].conj()

            write_buffer = A / (B - C - D)

            ## If no segment was folded in a given segment, make PSD a large float
            if len(missingSegments) > 0:
                write_buffer[:,missingSegments] = sys.float_info.max

            ifo1_psd             = np.sqrt(abs(write_buffer[:]))
            ifo2_psd             = np.sqrt(abs(write_buffer[:]))
              
            # Divide by the produce of adjact PSDs, so that when the pipeline divides
            # CSD by the product, it effectively processes data in which overlapping
            # window correction has already been acocunted for
            write_buffer_cplx = (write_buffer/winFactor) * foldedStatistic[iSegment:iSegment+nSegment_thisFile,:].conj()#.T

            #  If no segment was folded in a given segment, make CSD a small float
            if len(missingSegments) > 0:
                write_buffer_cplx[:,missingSegments] = sys.float_info.min
            if np.isreal(write_buffer_cplx).all():
                write_buffer_cplx = write_buffer_cplx.astype(complex)#np.complex(write_buffer_cplx)
                
            ReCSD                = np.real(write_buffer_cplx[:])
            ImCSD                = np.imag(write_buffer_cplx[:])

            del write_buffer
            del write_buffer_cplx
        
        #======================= Write to frames, at last! =======================#
        gpsStart = thisGPSStart
        segDuration = params['segmentDuration']
        siderealTime = siderealStart
        flow = params['flow']
        fhigh = params['fhigh']
        deltaF = params['deltaF']
        naiveBias = 1/((segDuration * deltaF * 2 - 1) * (9/11)) + 1
        bias = params['bias']
        w1w2bar = params['w1w2bar']
        w1w2squaredbar = params['w1w2squaredbar']
        w1w2ovlsquaredbar = params['w1w2ovlsquaredbar']
        nfreq = ((fhigh - flow) // deltaF) + 1
        datas = [{'name': ifo1+':LocalPSD', 'data': np.array(ifo1_psd), 'start': gpsStart, 'dx': foldedSegDuration/nfreq, 'kind': 'FSID'},
                 {'name': ifo2+':LocalPSD', 'data': np.array(ifo2_psd), 'start': gpsStart, 'dx': foldedSegDuration/nfreq, 'kind': 'FSID'},
                 {'name': ifo1+':AdjacentPSD', 'data': np.array(ifo1_psd), 'start': gpsStart, 'dx': foldedSegDuration/nfreq, 'kind': 'FSID'},
                 {'name': ifo2+':AdjacentPSD', 'data': np.array(ifo2_psd), 'start': gpsStart, 'dx': foldedSegDuration/nfreq, 'kind': 'FSID'},
                 {'name': ifo1+ifo2+':ReCSD', 'data': np.array(ReCSD), 'start': gpsStart, 'dx': foldedSegDuration/nfreq, 'kind': 'FSID'},
                 {'name': ifo1+ifo2+':ImCSD', 'data': np.array(ImCSD), 'start': gpsStart, 'dx': foldedSegDuration/nfreq, 'kind': 'FSID'},
                 {'name': ifo1+ifo2+':GPStime', 'data': np.array(gpsStart), 'start': gpsStart, 'dx': foldedSegDuration, 'kind': 'FSID'},
                 {'name': ifo1+ifo2+':SiderealTime', 'data': np.array(siderealTime), 'start': gpsStart, 'dx': foldedSegDuration, 'kind': 'FSID'},
                 {'name': ifo1+ifo2+':segmentDuration', 'data': np.array(segDuration), 'start': gpsStart, 'dx': foldedSegDuration, 'kind': 'FSID'},
                 {'name': ifo1+ifo2+':foldedSegmentDuration', 'data': np.array(foldedSegDuration), 'start': gpsStart, 'dx': foldedSegDuration, 'kind': 'FSID'},
                 {'name': ifo1+ifo2+':flow', 'data': np.array(flow), 'start': gpsStart, 'dx':foldedSegDuration, 'kind': 'FSID'},
                 {'name': ifo1+ifo2+':fhigh', 'data': np.array(fhigh), 'start': gpsStart, 'dx': foldedSegDuration, 'kind': 'FSID'},
                 {'name': ifo1+ifo2+':deltaF', 'data': np.array(deltaF), 'start': gpsStart, 'dx':foldedSegDuration, 'kind': 'FSID'},
                 {'name': ifo1+ifo2+':winFactor', 'data': np.array(winFactor), 'start': gpsStart, 'dx': foldedSegDuration, 'kind': 'FSID'},
                 {'name': ifo1+ifo2+':bias', 'data': np.array(bias), 'start': gpsStart, 'dx':foldedSegDuration, 'kind': 'FSID'},
                 {'name': ifo1+ifo2+':naivebias', 'data': np.array(naiveBias), 'start': gpsStart, 'dx': foldedSegDuration, 'kind': 'FSID'},
                 {'name': ifo1+ifo2+':w1w2bar', 'data': np.array(w1w2bar), 'start': gpsStart, 'dx': foldedSegDuration, 'kind': 'FSID'},
                 {'name': ifo1+ifo2+':w1w2squaredbar', 'data': np.array(w1w2squaredbar), 'start': gpsStart, 'dx':foldedSegDuration, 'kind': 'FSID'},
                 {'name': ifo1+ifo2+':w1w2ovlsquaredbar', 'data': np.array(w1w2ovlsquaredbar), 'start': gpsStart, 'dx':foldedSegDuration, 'kind': 'FSID'}]

        Fr.frputvect(frameFile,datas)

    # Last jobfile entry
    if (writeFSIDJobFile):
        if (jobStart >= 0):
            jobCount = jobCount + 1
            FSIDJobFilePtr.write(str(jobCount) + "\\" + str(jobStart) + '\\' + str(jobEnd) + '\\' + str(jobEnd-jobStart) + '\\')
            jobStart = -1 # Redundant
        FSIDJobFilePtr.close()
