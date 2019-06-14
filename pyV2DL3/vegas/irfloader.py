import ROOT
from root_numpy import hist2array
from ctypes import c_float
from pyV2DL3.vegas.load_vegas import VEGASStatus
from scipy.interpolate import RegularGridInterpolator
import numpy as np

def graph_to_array_y(graph):
    y = [g for g in graph.GetY()]
    return y


def graph_to_array_x(graph):
    x = [g for g in graph.GetX()]
    return x


def get_irf_not_safe(manager, offset_arr, az, ze, noise):
    offset_np_arr = np.array(offset_arr)
    loaded_offset = []
    ea_data_dict = {}
    ea = []
    ebias_data_dict = {}
    ebias = []
    abias_data_dict = {}
    abias = []
    for offset in offset_arr:
        effectiveAreaParameters = ROOT.VAEASimpleParameterData()
        effectiveAreaParameters.fAzimuth = az
        effectiveAreaParameters.fZenith = ze
        effectiveAreaParameters.fNoise = noise
        effectiveAreaParameters.fOffset = offset
        effectiveAreaParameters = manager.getVectorParamsFromSimpleParameterData(effectiveAreaParameters)
        ea_dl3 = manager.getEffectiveAreaCurve_DL3_no_theta_cut(effectiveAreaParameters)
        # Check if the effective area is empyt
        # NOTE: here we need to use == comparision
        # instead of the pythonic "is" operator 
        if(ea_dl3 == None):
            continue
        # Get Effective Area
        energy_bin = graph_to_array_x(ea_dl3)
        energy_bin = np.array(energy_bin)
        ea_array = graph_to_array_y(ea_dl3)
        ea_array = np.array(ea_array) 
        energyLow = np.power(10, energy_bin - (energy_bin[1] - energy_bin[0]) / 2.)
        energyHigh = np.power(10, energy_bin + (energy_bin[1] - energy_bin[0]) / 2.)
        ea.append(ea_array)
        ea_energy_low = energyLow
        ea_energy_high = energyHigh
        # Get Ebias
        a, e = hist2array(manager.getEnergyBias_DL3(effectiveAreaParameters, False),
                          return_edges=True)
        eLow = np.power(10, [e[0][:-1]])[0]
        eHigh = np.power(10, [e[0][1:]])[0]

        bLow = e[1][:-1]
        bHigh = e[1][1:]

        ac = []
        for aa in a:
            if np.sum(aa) > 0:
                ab = aa / (bHigh - bLow) / np.sum(aa)
            else:
                ab = aa
            try:
                ac = np.vstack((ac, ab))
            except:
                ac = ab

        ac = ac.transpose()
        ebias_energy_low = eLow
        ebias_energy_high = eHigh
        ebias_migration_low = bLow
        ebias_migration_high = bHigh
        ebias.append(ac)
        # Get ABias
        a, e = hist2array(manager.getAngularBias_DL3(effectiveAreaParameters),
                          return_edges=True)
        eLow = np.power(10, [e[0][:-1]])[0]
        eHigh = np.power(10, [e[0][1:]])[0]

        bLow = np.power(10, [e[1][:-1]])[0]
        bHigh = np.power(10, [e[1][1:]])[0]

        ac = []
        for aa in a:
            if np.sum(aa) > 0:
                # As the unit is sr^-1 we need to convert y bin size into radian  
                ab = aa / np.deg2rad(bHigh - bLow) / np.sum(aa) / np.pi / np.deg2rad(bHigh + bLow)
            else:
                ab = aa
            try:
                ac = np.vstack((ac, ab))
            except:
                ac = ab

        ac = ac.transpose()
        abias_energy_low = eLow
        abias_energy_high = eHigh
        abias_migration_low = bLow
        abias_migration_high = bHigh
        abias.append(ac)
        loaded_offset.append(offset)

    #  Create EA data
    ea = np.array(ea)
    ea_data_dict['ELow'] = ea_energy_low
    ea_data_dict['EHigh'] = ea_energy_high
    ea_data_dict['ThetaLow']  =  np.array(loaded_offset)#offset_np_arr
    ea_data_dict['ThetaHigh'] =  np.array(loaded_offset)#offset_np_arr
    ea_data_dict['Data'] = ea
    ebias = np.array(ebias)
    ebias_data_dict['ELow'] = ebias_energy_low
    ebias_data_dict['EHigh'] = ebias_energy_high
    ebias_data_dict['ThetaLow'] = np.array(loaded_offset)#offset_np_arr
    ebias_data_dict['ThetaHigh'] = np.array(loaded_offset)#offset_np_arr
    ebias_data_dict['MigrationLow'] = ebias_migration_low
    ebias_data_dict['MigrationHigh'] = ebias_migration_high
    ebias_data_dict['Data'] = ebias
    abias = np.array(abias)
    abias_data_dict['ELow'] = abias_energy_low
    abias_data_dict['EHigh'] = abias_energy_high
    abias_data_dict['ThetaLow'] = np.array(loaded_offset)#offset_np_arr
    abias_data_dict['ThetaHigh'] = np.array(loaded_offset)#offset_np_arr
    abias_data_dict['MigrationLow'] = abias_migration_low
    abias_data_dict['MigrationHigh'] = abias_migration_high
    abias_data_dict['Data'] = abias
    return ea_data_dict, ebias_data_dict, abias_data_dict


class IRFLoader:
    def __init__(self,vts_io):
        # If vegas not loaded. Load vegas
        self.__vegas__ = VEGASStatus()
        self.__vegas__.loadVEGAS()
        self.__manager__ = ROOT.VAEffectiveAreaManager()
        self.__manager__.loadEffectiveAreas(vts_io)
        self.__axis__ = ['Azimuth', 'Zenith', 'Noise', 'AbsoluteOffset']
        self.__buildIndex__()

    def __buildIndex__(self):
        manager = self.__manager__
        if (len(manager.fEffectiveAreas) <= 0):
            raise Exception("No effective areas! ")
        index_check = manager.fEffectiveAreas.at(0).fDimensionNames
        for k in self.__axis__:
            if (k not in index_check):
                raise Exception("IRF missing axis: {}".format(k))
        index_dict = {'Index': []}
        for i, ea in enumerate(manager.fEffectiveAreas):
            index_dict['Index'].append(i)
            for name, val in zip(ea.fDimensionNames, ea.fDimensionValues):
                if (not (name in index_dict.keys())):
                    index_dict[name] = []
                else:
                    index_dict[name].append(val)
        # Validate Completeness
        num_IRF = len(index_dict['Index'])
        axis_dict = {}
        check_num = 1
        for k in self.__axis__:
            check_num = check_num * len(np.unique(index_dict[k]))
            axis_dict[k] = np.sort(np.unique(index_dict[k]))
            if (len(axis_dict[k]) < 2):
                raise Exception('{} Axis need to have more than two values'.format(k))
        if (num_IRF != check_num):
            raise Exception("EA File not full")
        self.__axis_dict__ = axis_dict
        self.__index_dict__ = index_dict

    def getSafeEnergy(self,az,ze,noise):
        manager = self.__manager__
        effectiveAreaParameters = ROOT.VAEASimpleParameterData()
        effectiveAreaParameters.fAzimuth = az
        effectiveAreaParameters.fZenith = ze
        effectiveAreaParameters.fNoise = noise

        effectiveAreaParameters.fOffset = 0.5
        effectiveAreaParameters = manager.getVectorParamsFromSimpleParameterData(effectiveAreaParameters)
        minEnergy, maxEnergy = c_float(), c_float()
        # Is it the right way ? what does the offset here provide ?
        manager.getSafeEnergyRange(effectiveAreaParameters, 0.5, minEnergy, maxEnergy)
        return  minEnergy.value/1000.,maxEnergy.value/1000.

    def getIRF(self,az, ze, noise):
        # Find closest two values for az, ze and noise axis
        # Az
        az_index = self.__axis_dict__['Azimuth']

        if (az > 360):
            az_val_inter = az - 360
        else:
            az_val_inter = az
        for low, high in zip(az_index[:-1], az_index[1:]):
            if ((az >= low) and (az < high)):
                az_low = low
                az_high = high
                break
        if (az > az_index[-1]):
            az_low = az_index[-1]
            az_high = az_index[0] + 360
        # Ze
        ze_index = self.__axis_dict__['Zenith']
        ze_low = -1
        ze_high = -1
        for low, high in zip(ze_index[:-1], ze_index[1:]):
            if ((ze >= low) and (ze < high)):
                ze_low = low
                ze_high = high
                break
        if ((ze_low < 0) or (ze_high < 0)):
            raise Exception(" Ze out of range")
        # Noise
        noise_index = self.__axis_dict__['Noise']
        noise_low = -1
        noise_high = -1
        for low, high in zip(noise_index[:-1], noise_index[1:]):
            if ((noise >= low) and (noise < high)):
                noise_low = low
                noise_high = high
                break
        if ((noise_low < 0) or (noise_high < 0)):
            raise Exception("Noise out of range")
        # Done finding index values use for interpolation

        irf_data = []
        offset_index = self.__axis_dict__['AbsoluteOffset']
        for az_i, az_val in [(0, az_low), (1, az_high)]:
            for ze_i, ze_val in[(0, ze_low), (1, ze_high)]:
                for noise_i, noise_val in [(0, noise_low), (1, noise_high)]:
                    irf_dict = {}
                    irf_dict['Index'] = (az_i, ze_i, noise_i)
                    az_shifted = az_val if az_val <= 360 else az_val - 360.
                    ea_dict, ebias_dict, abias_dict = get_irf_not_safe(self.__manager__,
                                                                       offset_index,
                                                                       az_shifted,
                                                                       ze_val,
                                                                       noise_val)
                    irf_dict['EA_Dict'] = ea_dict
                    irf_dict['EBias_Dict'] = ebias_dict
                    irf_dict['ABias_Dict'] = abias_dict
                    irf_data.append(irf_dict)
                # Load values
        # Initilaize Data container
        ea_edim = 0
        elow = [] 
        offset_dim = 0
        offset_low = []
        for irf in irf_data:
            ea_data_peek = irf['EA_Dict']['Data']
            ebias_data_peek = irf['EBias_Dict']['Data']
            abias_data_peek = irf['ABias_Dict']['Data']
            if(ea_edim < ea_data_peek.shape[1]):
                ea_edim = ea_data_peek.shape[1]
                elow   = irf['EA_Dict']['ELow'] 
            if(offset_dim < ea_data_peek.shape[0]):
                offset_dim = ea_data_peek.shape[0]
                offset_low = irf['EA_Dict']['ThetaLow']
        
        ea_array = np.zeros([2, 2, 2,offset_dim,ea_edim ])
        ebias_array = np.zeros([2, 2, 2,offset_dim, (ebias_data_peek.shape[1]),(ebias_data_peek.shape[1])])
        abias_array = np.zeros([2, 2, 2,offset_dim, (abias_data_peek.shape[1]),(abias_data_peek.shape[1])])

        # Build Interpolator
        for irf in irf_data:
            index = irf['Index']
            ea_data = irf['EA_Dict']['Data']
            ebias_data = irf['EBias_Dict']['Data']
            abias_data = irf['ABias_Dict']['Data']
            elow_first = irf['EA_Dict']['ELow'][0] 
            offset_low_first = irf['EA_Dict']['ThetaLow'][0]

            for i,val in enumerate(elow):
              if(np.abs(val -  elow_first) < 1e-15): 
                break 
            eindex_low,eindex_high = i,i+len(irf['EA_Dict']['ELow'])

            for i,val in enumerate(offset_low):
               if(val == offset_low_first):
                break 
            offset_index_low,offset_index_high = i,i+len(irf['EA_Dict']['ThetaLow'])

            ea_array[index[0], index[1], index[2]][offset_index_low:offset_index_high,eindex_low:eindex_high] = ea_data
            ebias_array[index[0], index[1], index[2]][offset_index_low:offset_index_high] = ebias_data
            abias_array[index[0], index[1], index[2]][offset_index_low:offset_index_high]= abias_data
        inter_axis = np.array([[az_low, az_high],
                               [ze_low, ze_high],
                               [noise_low, noise_high]])
        ea_interpolator = RegularGridInterpolator(inter_axis, ea_array)
        ebias_interpolator = RegularGridInterpolator(inter_axis, ebias_array)
        abias_interpolator = RegularGridInterpolator(inter_axis, abias_array)

        # Now lets actually build the data block to be passed
        # EA
        elow = irf_data[0]['EA_Dict']['ELow']
        ehigh = irf_data[0]['EA_Dict']['EHigh']
        thetalow = irf_data[0]['EA_Dict']['ThetaLow']
        thetahigh = irf_data[0]['EA_Dict']['ThetaHigh']
        ea_interpolated = ea_interpolator((az, ze, noise))
        ea_final_data = np.array([(elow, ehigh, thetalow, thetahigh, ea_interpolated)],
                                 dtype=[('ENERG_LO', '>f4', np.shape(elow)),
                                        ('ENERG_HI', '>f4', np.shape(ehigh)),
                                        ('THETA_LO', '>f4', np.shape(thetalow)),
                                        ('THETA_HI', '>f4', np.shape(thetahigh)),
                                        ('EFFAREA', '>f4', np.shape(ea_interpolated))])
        # EBias
        elow = irf_data[0]['EBias_Dict']['ELow']
        ehigh = irf_data[0]['EBias_Dict']['EHigh']
        thetalow = irf_data[0]['EBias_Dict']['ThetaLow']
        thetahigh = irf_data[0]['EBias_Dict']['ThetaHigh']
        miglow = irf_data[0]['EBias_Dict']['MigrationLow']
        mighigh = irf_data[0]['EBias_Dict']['MigrationHigh']
        ebias_interpolated = ebias_interpolator((az, ze, noise))

        ebias_final_data = np.array([(elow, ehigh, miglow, mighigh, thetalow, thetahigh, ebias_interpolated)],
                                    dtype=[('ENERG_LO', '>f4', np.shape(elow)),
                                           ('ENERG_HI', '>f4', np.shape(ehigh)),
                                           ('MIGRA_LO', '>f4', np.shape(miglow)),
                                           ('MIGRA_HI', '>f4', np.shape(mighigh)),
                                           ('THETA_LO', '>f4', np.shape(thetalow)),
                                           ('THETA_HI', '>f4', np.shape(thetahigh)),
                                           ('MATRIX', '>f4', np.shape(ebias_interpolated))])
        # ABias
        elow = irf_data[0]['ABias_Dict']['ELow']
        ehigh = irf_data[0]['ABias_Dict']['EHigh']
        thetalow = irf_data[0]['ABias_Dict']['ThetaLow']
        thetahigh = irf_data[0]['ABias_Dict']['ThetaHigh']
        miglow   = irf_data[0]['ABias_Dict']['MigrationLow']
        mighigh   = irf_data[0]['ABias_Dict']['MigrationHigh']
        abias_interpolated = abias_interpolator((az, ze, noise))
        # Flip axis order
        # Axis order:
        # Energy, Theta, Rad
        abias_interpolated = np.transpose(abias_interpolated, axes=(1, 0, 2))
        abias_final_data = np.array([(elow, ehigh, thetalow, thetahigh, miglow, mighigh, abias_interpolated)],
                                    dtype=[('ENERG_LO', '>f4', np.shape(elow)),
                                           ('ENERG_HI', '>f4', np.shape(ehigh)),
                                           ('THETA_LO', '>f4', np.shape(thetalow)),
                                           ('THETA_HI', '>f4', np.shape(thetahigh)),
                                           ('RAD_LO', '>f4', np.shape(miglow)),
                                           ('RAD_HI', '>f4', np.shape(mighigh)),
                                           ('RPSF', '>f4', np.shape(abias_interpolated))])
        return ea_final_data, ebias_final_data, abias_final_data