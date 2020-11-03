import uproot
import numpy as np
import logging
from pyV2DL3.eventdisplay.util import bin_centers_to_edges
from pyV2DL3.eventdisplay.IrfInterpolator import IrfInterpolator

logger = logging.getLogger(__name__)


def __fillRESPONSE_not_safe__(effectiveArea, azimuth, zenith, noise, offset, irf_to_store={}):
    response_dict = {}
    from ROOT import VPlotInstrumentResponseFunction, VGammaHadronCuts, TTree

    filename = effectiveArea.GetName()
    cuts = effectiveArea.Get("GammaHadronCuts")
    # EventDisplay IRF interpolator object
    irf_interpolator = IrfInterpolator(filename, azimuth)

    # Extract the camera offsets simulated within the effective areas file.
    fast_eff_area = uproot.open(filename)['fEffArea']
    camera_offsets = np.unique(np.round(fast_eff_area.array('Woff'), decimals=2))
    offset = camera_offsets  # new: to add all offsets
    print("camera offsets: ", camera_offsets, offset, "noise:",noise,"zenith:",zenith)
    # Check the camera offset bins available in the effective area file.
    theta_low = []
    theta_high = []
    if len(camera_offsets) == 1:
        # Many times, just IRFs for 0.5 deg. Assume that offset for the whole camera.
        logger.debug('IMPORTANT: Only one camera offset bin '
                     + '({} deg) simulated within the effective area file selected.'.format(camera_offsets[0]))
        logger.debug('IMPORTANT: Setting the IRFs of that given camera offset value to the whole camera')
        theta_low = [0.0, 10.0]
        theta_high = [0.0, 10.0]
    if len(camera_offsets) > 1:
        # Note in the camera offset _low and _high may refer to the simulated "points", and
        # not to actual bins.
        theta_low = camera_offsets
        theta_high = camera_offsets
        print("theta low inside if clause if camera offset >1", theta_low)
    print(theta_low)
    logger.debug('Getting Theta2 cut from EA file')
    theta2cut = cuts.fCut_Theta2_max
    logger.debug('Theta2 cut is {:.2f}'.format(theta2cut))

    if irf_to_store['point-like']:
        #
        # Interpolate effective area  (point-like)
        #
        irf_interpolator.set_irf('eff')
        eff_area, axis = irf_interpolator.interpolate([noise, zenith, offset])
        log_energy_tev = axis[0]
        energy_low = np.power(10, log_energy_tev - (log_energy_tev[1] - log_energy_tev[0]) / 2.)
        energy_high = np.power(10, log_energy_tev + (log_energy_tev[1] - log_energy_tev[0]) / 2.)

        y = np.array(eff_area)
        ea = [y, y]

        x = np.array([(energy_low, energy_high, theta_low, theta_high, ea)],
                     dtype=[('ENERG_LO', '>f4', np.shape(energy_low)),
                            ('ENERG_HI', '>f4', np.shape(energy_high)),
                            ('THETA_LO', '>f4', np.shape(theta_low)),
                            ('THETA_HI', '>f4', np.shape(theta_high)),
                            ('EFFAREA', '>f4', np.shape(ea))])
        response_dict['EA'] = x
        response_dict['LO_THRES'] = min(energy_low)
        response_dict['HI_THRES'] = max(energy_high)
        response_dict['RAD_MAX'] = np.sqrt(theta2cut)
        #
        # Energy dispersion (point-like)
        #
        irf_interpolator.set_irf('hEsysMCRelative2D')
        bias, axis = irf_interpolator.interpolate([noise, zenith, offset])

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
        print("ac:", ac)
        print(len(theta_low))
        x = np.array([(eLow, eHigh, bLow, bHigh, theta_low, theta_high, [ac, ac])],
                     dtype=[('ENERG_LO', '>f4', (len(eLow),)),
                            ('ENERG_HI', '>f4', (len(eHigh),)),
                            ('MIGRA_LO', '>f4', (len(bLow),)),
                            ('MIGRA_HI', '>f4', (len(bHigh),)),
                            ('THETA_LO', '>f4', (len(theta_low),)),
                            ('THETA_HI', '>f4', (len(theta_high),)),
                            ('MATRIX', '>f4', (len(theta_low), np.shape(ac)[0], np.shape(ac)[1]))])
        print('before')
        print(x)
        response_dict['MIGRATION'] = x
        response_dict['RAD_MAX'] = np.sqrt(theta2cut)
        print('IRF interpolation done')
    if irf_to_store['full-enclosure']:
        #
        # Interpolate effective area (full-enclosure)
        #
        # irf_interpolator.set_irf('gEffAreaNoTh2MC')
        irf_interpolator.set_irf('effNoTh2')
        eff_area, axis = irf_interpolator.interpolate([noise, zenith, offset])
        log_energy_tev = axis[0]
        energy_low = np.power(10, log_energy_tev - (log_energy_tev[1] - log_energy_tev[0]) / 2.)
        energy_high = np.power(10, log_energy_tev + (log_energy_tev[1] - log_energy_tev[0]) / 2.)

        y = np.array(eff_area)
        ea = [y, y]

        x = np.array([(energy_low, energy_high, theta_low, theta_high, ea)],
                     dtype=[('ENERG_LO', '>f4', np.shape(energy_low)),
                            ('ENERG_HI', '>f4', np.shape(energy_high)),
                            ('THETA_LO', '>f4', np.shape(theta_low)),
                            ('THETA_HI', '>f4', np.shape(theta_high)),
                            ('EFFAREA', '>f4', np.shape(ea))])
        response_dict['LO_THRES'] = min(energy_low)
        response_dict['HI_THRES'] = max(energy_high)
        response_dict['FULL_EA'] = x
        #
        # Energy dispersion (full-enclosure)
        #
        irf_interpolator.set_irf('hEsysMCRelative2DNoDirectionCut')
        print(noise, zenith, offset)
        bias, axis = irf_interpolator.interpolate([noise, zenith, offset])

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
        print ("len, :eLow, eHigh, bLow, bHigh, theta_low, theta_high,ac",(len(eLow),),(len(eHigh),),(len(bLow),),(len(bHigh),),(len(theta_low),),(len(theta_high),),np.shape(ac)[0],np.shape(ac)[1])
        # print ([ac,ac])
        # print ("len [ac,ac]",len([ac,ac],))
        x = np.array([(eLow, eHigh, bLow, bHigh, theta_low, theta_high, [ac, ac])],
                     dtype=[('ENERG_LO', '>f4', (len(eLow),)),
                            ('ENERG_HI', '>f4', (len(eHigh),)),
                            ('MIGRA_LO', '>f4', (len(bLow),)),
                            ('MIGRA_HI', '>f4', (len(bHigh),)),
                            ('THETA_LO', '>f4', (len(theta_low),)),
                            ('THETA_HI', '>f4', (len(theta_high),)),
                            ('MATRIX', '>f4', (len(theta_low), np.shape(ac)[0], np.shape(ac)[1]))])

        response_dict['FULL_MIGRATION'] = x
        #
        # Direction dispersion (for full-enclosure IRFs)
        #
        irf_interpolator.set_irf('hAngularLogDiffEmc_2D')
        direction_diff, axis = irf_interpolator.interpolate([noise, zenith, offset])

        energy_edges = bin_centers_to_edges(axis[0])
        rad_edges = bin_centers_to_edges(axis[1])

        eLow = np.power(10, [energy_edges[:-1]])[0]
        eHigh = np.power(10, [energy_edges[1:]])[0]

        rLow = np.power(10, [rad_edges[:-1]])[0]
        rHigh = np.power(10, [rad_edges[1:]])[0]

        ac = []
        for aa in direction_diff.transpose():
            if np.sum(aa) > 0:
                ab = aa / np.sum(aa * (rHigh - rLow))
            else:
                ab = aa
            try:
                ac = np.vstack((ac, ab))
            except:
                ac = ab

        ac = ac.transpose()
        x = np.array([(eLow, eHigh, rLow, rHigh, theta_low, theta_high, [ac, ac])],
                     dtype=[('ENERG_LO', '>f4', (len(eLow),)),
                            ('ENERG_HI', '>f4', (len(eHigh),)),
                            ('RAD_LO', '>f4', (len(rLow),)),
                            ('RAD_HI', '>f4', (len(rHigh),)),
                            ('THETA_LO', '>f4', (len(theta_low),)),
                            ('THETA_HI', '>f4', (len(theta_high),)),
                            ('RPSF', '>f4', (len(theta_low), np.shape(ac)[0], np.shape(ac)[1]))])

        response_dict['PSF'] = x
    return response_dict
