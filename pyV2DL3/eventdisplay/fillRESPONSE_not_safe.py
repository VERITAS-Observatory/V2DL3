from astropy.io import fits
from ctypes import c_float
import numpy as np
import logging
from pyV2DL3.eventdisplay.util import is_close
from root_numpy import hist2array
import ROOT
 
logger = logging.getLogger(__name__)


def __fillRESPONSE_not_safe__(effectiveArea, azimuth, zenith, noise, offset):
    response_dict = {}

    from ROOT import VPlotInstrumentResponseFunction, GammaHadronCuts, TTree

    filename = effectiveArea.GetName()
    cuts = effectiveArea.Get("GammaHadronCuts")
    effAreaTree = effectiveArea.Get("fEffArea")

    # FIXME: It's sad, but this is the only way I find to extract the IRFs...
    plotter = VPlotInstrumentResponseFunction()
    plotter.addInstrumentResponseData(filename, 0., 0.5, 0, 1.5, 200, "A_MC")
    # Get Theta2 cut from file
    logger.debug('Getting Theta2 cut from EA file')
    theta2cut = cuts.fCut_Theta2_max
    logger.debug('Theta2 cut is {:.2f}'.format(theta2cut))

    # Extract effective area:
    # FIXME: For now, just manually extracting the correct value for the 3ML runs
    # avAlt = 80.17160359513885
    # avAz = -170.57748918711198
    # index = 2.3
    # pedVars = 5.438850227591036         # Closest value: 5.22773 noise = 200
    # FIXME: We choose azimuth bin between [-180.0, -120.0]
    for entry in effAreaTree:
        # Now you have acess to the leaves/branches of each entry in the tree, e.g.
        #     if (np.abs(closest_ze - entry.ze) < diff_ze):
        #         closest_ze = entry.ze
        if is_close(entry.ze, 20.) and is_close(entry.azMin, -180.) and is_close(entry.azMax, -120.) and is_close(
                entry.index, 2.3):
            effArea = np.array([i for i in entry.eff])
            log_energy_TeV = np.array([i for i in entry.e0])
            bias, en = hist2array(entry.hEsysMCRelative2D, return_edges=True)
            break
            # print(entry.ze, entry.azMin, entry.azMax, entry.Woff, entry.noise, entry.pedvar, entry.index)
    energyLow = np.power(10, log_energy_TeV - (log_energy_TeV[1] - log_energy_TeV[0]) / 2.)
    energyHigh = np.power(10, log_energy_TeV + (log_energy_TeV[1] - log_energy_TeV[0]) / 2.)
    thetaLow = [0.0, 1.0]
    thetaHigh = [1.0, 2.0]


    # effectiveAreaManager = ROOT.VAEffectiveAreaManager()
    # effectiveAreaManager.loadEffectiveAreas(effectiveAreaIO)
    # effectiveAreaManager.setUseReconstructedEnergy(False)
    #
    # effectiveAreaParameters = ROOT.VAEASimpleParameterData()
    # effectiveAreaParameters.fAzimuth = (azimuth)
    # effectiveAreaParameters.fZenith = zenith
    # effectiveAreaParameters.fNoise = noise
    # effectiveAreaParameters.fOffset = offset
    # effectiveArea = effectiveAreaManager.getEffectiveAreaCurve(effectiveAreaParameters)
    # effectiveAreaParameters = effectiveAreaManager.getVectorParamsFromSimpleParameterData(effectiveAreaParameters)
    # # Filling effective area
    # x, y, ye = [], [], []
    # for i in range(effectiveArea.GetN()):
    #     tmpX, tmpY = ROOT.Double(0), ROOT.Double(0)
    #     effectiveArea.GetPoint(i, tmpX, tmpY)
    #     ye.append(effectiveArea.GetErrorY(i))
    #     effectiveArea.GetErrorX
    #     x.append(tmpX)
    #     y.append(tmpY)
    #
    # x = np.array(x)
    # y = np.array(y)
    # ye = np.array(ye)
    # energyLow = np.power(10, x - (x[1] - x[0])/2.)
    # energyHigh = np.power(10, x + (x[1] - x[0])/2.)
    # thetaLow = [0.0, 1.0]
    # thetaHigh = [1.0, 2.0]
    # ea = np.vstack((y, y))
    y = np.array(effArea)
    ea = [y, y]
    # minEnergy , maxEnergy = c_float(), c_float()
    # effectiveAreaManager.getSafeEnergyRange(effectiveAreaParameters, 0.5, minEnergy, maxEnergy)

    x = np.array([(energyLow, energyHigh, thetaLow, thetaHigh, ea)], 
                 dtype=[('ENERG_LO', '>f4', np.shape(energyLow)), 
                        ('ENERG_HI', '>f4', np.shape(energyHigh)), 
                        ('THETA_LO', '>f4', np.shape(thetaLow)), 
                        ('THETA_HI', '>f4', np.shape(thetaHigh)), 
                        ('EFFAREA', '>f4', np.shape(ea))])
    response_dict['EA'] = x
    response_dict['LO_THRES'] = min(energyLow)
    response_dict['HI_THRES'] = max(energyHigh)
    response_dict['RAD_MAX'] = np.sqrt(theta2cut)

    eLow = np.power(10, [en[0][:-1]])[0]
    eHigh = np.power(10, [en[0][1:]])[0]

    bLow = np.array([en[1][:-1]])[0]
    bHigh = np.array([en[1][1:]])[0]
    
    ac = []
    for aa in bias:
        if np.sum(aa) > 0:
            ab = aa / np.sum(aa*(bHigh - bLow))
        else:
            ab = aa
        try:
            ac = np.vstack((ac, ab))
        except:
            ac = ab
            
    ac = ac.transpose() 
    x = np.array([(eLow, eHigh, bLow, bHigh, [0, 1.0], [1.0, 2.0], [ac, ac])], 
                 dtype=[('ETRUE_LO', '>f4', (len(eLow),)), 
                        ('ETRUE_HI', '>f4', (len(eHigh),)), 
                        ('MIGRA_LO', '>f4', (len(bLow),)), 
                        ('MIGRA_HI', '>f4', (len(bLow),)), 
                        ('THETA_LO', '>f4', (2,)), 
                        ('THETA_HI', '>f4', (2,)), 
                        ('MATRIX', '>f4', (2, np.shape(ac)[0], np.shape(ac)[1]))])
    
    response_dict['MIGRATION'] = x
    return response_dict 
