from astropy.io import fits
from ctypes import c_float
import numpy as np
import logging
from pyV2DL3.eventdisplay.util import bin_centers_to_edges
from root_numpy import hist2array
from pyV2DL3.eventdisplay.IrfInterpolator import IrfInterpolator
import ROOT

logger = logging.getLogger(__name__)


def __fillRESPONSE_not_safe__(effectiveArea, azimuth, zenith, noise, offset):
    response_dict = {}

    from ROOT import VPlotInstrumentResponseFunction, GammaHadronCuts, TTree

    filename = effectiveArea.GetName()
    cuts = effectiveArea.Get("GammaHadronCuts")
    effAreaTree = effectiveArea.Get("fEffArea")
    # EventDisplay IRF interpolator object
    irf_interpolator = IrfInterpolator(filename)

    # FIXME: This whole part of the code needs to be updated with the new 'IrfInterpolator' class.
    # plotter = VPlotInstrumentResponseFunction()
    # plotter.addInstrumentResponseData(filename, 0., 0.5, 0, 1.5, 200, "A_MC")
    # Get Theta2 cut from file
    logger.debug('Getting Theta2 cut from EA file')
    theta2cut = cuts.fCut_Theta2_max
    logger.debug('Theta2 cut is {:.2f}'.format(theta2cut))
    #
    # Interpolate effective area:
    #
    irf_interpolator.set_irf('eff')
    eff_area, axis = irf_interpolator.interpolate([0, 3., noise, zenith, offset])
    log_energy_TeV = axis[0]
    energyLow = np.power(10, log_energy_TeV - (log_energy_TeV[1] - log_energy_TeV[0]) / 2.)
    energyHigh = np.power(10, log_energy_TeV + (log_energy_TeV[1] - log_energy_TeV[0]) / 2.)

    # Extract effective area:
    # FIXME: Choose the best azimuth bin
    # FIXME: Check the camera offset, and add all available bins
    thetaLow = [0.0, 10.0]
    thetaHigh = [0.0, 10.0]

    y = np.array(eff_area)
    ea = [y, y]

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
    #
    # Energy dispersion
    #
    irf_interpolator.set_irf('hEsysMCRelative2D')
    bias, axis = irf_interpolator.interpolate([0, 3., noise, zenith, offset])

    energy_edges = bin_centers_to_edges(axis[0])
    bias_edges = bin_centers_to_edges(axis[1])

    eLow = np.power(10, [energy_edges[:-1]])[0]
    eHigh = np.power(10, [energy_edges[1:]])[0]

    bLow = np.array([bias_edges[:-1]])[0]
    bHigh = np.array([bias_edges[1:]])[0]

    ac = []
    for aa in bias.transpose():
        if np.sum(aa) > 0:
            ab = aa / np.sum(aa * (bHigh - bLow))
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
