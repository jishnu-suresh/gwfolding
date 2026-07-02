def loadSIDParams(SIDFile, ifo1, ifo2, segDuration):
    #
    # Load SID parameters from a frame
    #
    #
    # param = loadSIDParams (SIDFile, ifo1, ifo2);
    #
    # SIDFile      : String: File containing the SID
    # ifo1, ifo2   : Strings: Detectors in the baseline
    # segDuration  : Real/Integer: Suration of each segment (possibly 52sec)
    #
    # param        : Dictionary: Different output parameters
    #

    # 
    # Author: Sanjit Mitra <sanjit.mitra@ligo.org>
    # Translated to Python by Erik Floden <erik.floden@ligo.org>
    from Fr import frgetvect as gv

    # Insert parameters in the output structure
    param = {}
    param['segmentDuration'] = gv(SIDFile,ifo1+ifo2+':segmentDuration',0,segDuration)[0][0]
    param['flow'] = gv(SIDFile,ifo1+ifo2+':flow',0,segDuration)[0][0]
    param['fhigh'] = gv(SIDFile,ifo1+ifo2+':fhigh',0,segDuration)[0][0]
    param['deltaF'] = gv(SIDFile,ifo1+ifo2+':deltaF',0,segDuration)[0][0]
    param['w1w2bar'] = gv(SIDFile,ifo1+ifo2+':w1w2bar',0,segDuration)[0][0]
    param['w1w2squaredbar'] = gv(SIDFile,ifo1+ifo2+':w1w2squaredbar',0,segDuration)[0][0]
    param['w1w2ovlsquaredbar'] = gv(SIDFile,ifo1+ifo2+':w1w2ovlsquaredbar',0,segDuration)[0][0]
    param['bias'] = gv(SIDFile,ifo1+ifo2+':bias',0,segDuration)[0][0]
    
    return param