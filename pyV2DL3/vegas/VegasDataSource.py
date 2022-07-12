import ROOT

from pyV2DL3.vegas.fillEVENTS_not_safe import __fillEVENTS_not_safe__
from pyV2DL3.vegas.fillRESPONSE_not_safe import __fillRESPONSE_not_safe__
from pyV2DL3.vegas.load_vegas import VEGASStatus
from pyV2DL3.VtsDataSource import VtsDataSource


class VegasDataSource(VtsDataSource):
    def __init__(self, etv_file, ea_file,
                 event_classes=None,
                 save_msw_msl=False,
                 ):
        super(VegasDataSource, self).__init__("VEGAS", etv_file, ea_file)

        # Developer exceptions to ensure this was constructed with an EA xor EventClass(es)
        if ea_file is None and event_classes is None:
            raise Exception("Running V2DL3 without effective area file(s) is currently unsupported.")
        elif ea_file is not None and event_classes is not None:
            raise Exception("VegasDataSource was somehow constructed with both an effective area"
                            + " file and event class")

        # Loading VEGAS if not already done so
        self.vegas_status = VEGASStatus()
        self.vegas_status.loadVEGAS()
        self.__evt_file__ = ROOT.VARootIO(etv_file, True)
        self.__event_classes__ = event_classes
        self.__save_msw_msl__ = save_msw_msl

        if ea_file is not None:
            self.__ea_file__ = ROOT.VARootIO(ea_file, True)

        # Auxiliary storage
        self.__azimuth__ = 0
        self.__zenith__ = 0
        self.__noise__ = 0

    def __del__(self):
        """Close the root files

        These typechecks will prevent the user from having their true exception
        buried by a CPyCppyy exception on program exit.
        """
        cpy_nonestring = "<class 'CPyCppyy_NoneType'>"

        if str(type(self.__evt_file__)) != cpy_nonestring and not isinstance(self.__evt_file__, str):
            self.__evt_file__.closeTheRootFile()

        if self.__ea_file__ is not None:
            if str(type(self.__ea_file__)) != cpy_nonestring and not isinstance(self.__ea_file__, str):
                self.__ea_file__.closeTheRootFile()

    def __fill_evt__(self):
        gti, ea_config, evt_dicts = __fillEVENTS_not_safe__(self.__evt_file__,
                                                            event_classes=self.__event_classes__,
                                                            save_msw_msl=self.__save_msw_msl__,
                                                            )
        self.__gti__ = gti
        # This is an array of dicts for each event class (array of one when not using event class mode)
        self.__evt__ = evt_dicts
        self.__azimuth__ = ea_config["azimuth"]
        self.__zenith__ = ea_config["zenith"]
        self.__noise__ = ea_config["noise"]

    def __fill_gti__(self):
        pass

    def __fill_response__(self):
        az = self.__azimuth__
        ze = self.__zenith__
        nn = self.__noise__
        response_dicts = []
        # Event class mode
        if self.__event_classes__ is not None:
            # Fill responses per event class
            for ec in self.__event_classes__:
                ea = ec.effective_area_IO
                response_dicts.append(
                    __fillRESPONSE_not_safe__(ea, az, ze, nn, self.__irf_to_store__)
                )
        # When not using event class mode
        else:
            response_dicts.append(
                __fillRESPONSE_not_safe__(self.__ea_file__, az, ze, nn, self.__irf_to_store__)
            )
        self.__response__ = response_dicts
