from astropy.io import fits
from ctypes import c_float
import numpy as np
import logging
from root_numpy import hist2array
from pyV2DL3.vegas.util import decodeConfigMask,getThetaSquareCut
import ROOT
 
logger = logging.getLogger(__name__)

def __fillRESPONSE_not_safe__(effectiveAreaIO,azimuth,zenith,noise,offset):
    response_dict = {}

    effectiveAreaIO.loadTheRootFile() 
    # Get Theta2 cut from file
    logger.debug('Getting Theta2 cut from EA file')
    cuts = effectiveAreaIO.loadTheCutsInfo()
    for k in cuts:
        theta2cut = getThetaSquareCut(k.fCutsFileText)

    logger.debug('Theta2 cut is {:.2f}'.format(theta2cut))
    effectiveAreaManager = ROOT.VAEffectiveAreaManager()
    effectiveAreaManager.loadEffectiveAreas(effectiveAreaIO)
    effectiveAreaManager.setUseReconstructedEnergy(False)    

    effectiveAreaParameters = ROOT.VAEASimpleParameterData()
    effectiveAreaParameters.fAzimuth = (azimuth)
    effectiveAreaParameters.fZenith = zenith 
    effectiveAreaParameters.fNoise = noise 
    effectiveAreaParameters.fOffset = offset 
    effectiveArea = effectiveAreaManager.getEffectiveAreaCurve(effectiveAreaParameters)
    effectiveAreaParameters = effectiveAreaManager.getVectorParamsFromSimpleParameterData(effectiveAreaParameters)
    # Filling effective area
    x, y, ye = [], [], []
    for i in range(effectiveArea.GetN()):
        tmpX, tmpY = ROOT.Double(0), ROOT.Double(0)
        effectiveArea.GetPoint(i, tmpX, tmpY)
        ye.append(effectiveArea.GetErrorY(i))
        effectiveArea.GetErrorX
        x.append(tmpX)
        y.append(tmpY)
            
    x = np.array(x)
    y = np.array(y)
    ye = np.array(ye)
    energyLow = np.power(10, x - (x[1] - x[0])/2.)
    energyHigh = np.power(10, x + (x[1] - x[0])/2.)
    thetaLow = [0.0, 1.0]
    thetaHigh = [1.0, 2.0]
    # ea = np.vstack((y, y))
    ea = [y,y]
    minEnergy , maxEnergy = c_float(), c_float()
    effectiveAreaManager.getSafeEnergyRange(effectiveAreaParameters, 0.5, minEnergy, maxEnergy)

    x = np.array([(energyLow, energyHigh, thetaLow, thetaHigh, ea)], 
                 dtype=[('ENERG_LO', '>f4', np.shape(energyLow)), 
                        ('ENERG_HI', '>f4', np.shape(energyHigh)), 
                        ('THETA_LO', '>f4', np.shape(thetaLow)), 
                        ('THETA_HI', '>f4', np.shape(thetaHigh)), 
                        ('EFFAREA', '>f4', np.shape(ea))])
    response_dict['EA'] = x
    response_dict['LO_THRES'] = minEnergy.value/1000.
    response_dict['HI_THRES'] = maxEnergy.value/1000.
    response_dict['RAD_MAX']  = np.sqrt(theta2cut)

    a, e = hist2array(effectiveAreaManager.getEnergyBias2D(effectiveAreaParameters), return_edges=True)
    eLow = np.power(10, [e[0][:-1]])[0]
    eHigh = np.power(10, [e[0][1:]])[0]
    
    bLow = np.power(10, [e[1][:-1]])[0]
    bHigh = np.power(10, [e[1][1:]])[0]
    
    ac = []
    for aa in a:
        if np.sum(aa) > 0:
            ab = aa / np.sum(aa*(bHigh - bLow))
        else:
            ab = aa
        try:
            ac = np.vstack((ac, ab))
        except:
            ac = ab
            
    ac = ac.transpose() 
    x = np.array([(eLow, eHigh, bLow, bHigh, [0, 10.], [0, 10.], [ac, ac])],
                 dtype=[('ENERG_LO', '>f4', (len(eLow),)), 
                        ('ENERG_HI', '>f4', (len(eHigh),)), 
                        ('MIGRA_LO', '>f4', (len(bLow),)), 
                        ('MIGRA_HI', '>f4', (len(bLow),)), 
                        ('THETA_LO', '>f4', (2,)), 
                        ('THETA_HI', '>f4', (2,)), 
                        ('MATRIX', '>f4', (2, np.shape(ac)[0], np.shape(ac)[1]))])
    
    response_dict['MIGRATION'] = x
    return response_dict 
